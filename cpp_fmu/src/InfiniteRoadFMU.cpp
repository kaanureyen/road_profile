#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "fmi2Functions.h"

// Custom random generator matching numpy.random.RandomState(seed)
class NumPyRandom {
private:
    unsigned long mt[624];
    int mti;

public:
    NumPyRandom(unsigned long seed) {
        mt[0] = seed & 0xffffffffUL;
        for (mti = 1; mti < 624; mti++) {
            mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            mt[mti] &= 0xffffffffUL;
        }
    }

    unsigned long next() {
        unsigned long y;
        static unsigned long mag01[2] = {0x0UL, 0x9908b0dfUL};

        if (mti >= 624) {
            int kk;
            for (kk = 0; kk < 624 - 397; kk++) {
                y = (mt[kk] & 0x80000000UL) | (mt[kk + 1] & 0x7fffffffUL);
                mt[kk] = mt[kk + 397] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            for (; kk < 624 - 1; kk++) {
                y = (mt[kk] & 0x80000000UL) | (mt[kk + 1] & 0x7fffffffUL);
                mt[kk] = mt[kk + (397 - 624)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            y = (mt[624 - 1] & 0x80000000UL) | (mt[0] & 0x7fffffffUL);
            mt[624 - 1] = mt[397 - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mti = 0;
        }

        y = mt[mti++];

        /* Tempering */
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        return y;
    }

    double uniform(double low, double high) {
        unsigned long a = next() >> 5;
        unsigned long b = next() >> 6;
        double r = (a * 67108864.0 + b) / 9007199254740992.0;
        return low + r * (high - low);
    }
};

struct ModelInstance {
    std::string instanceName;
    std::string guid;
    fmi2CallbackLogger logger = nullptr;
    fmi2ComponentEnvironment componentEnvironment = nullptr;
    fmi2Boolean loggingOn = fmi2False;

    // Inputs (Real)
    double x = 0.0;
    double y = 0.0;

    // Outputs (Real)
    double z = 0.0;

    // Parameters
    int seed = 42;
    int road_class = 3; // 1=A, 2=B, 3=C, 4=D, 5=E, 0=Custom Gd_n0
    double Gd_n0 = 256e-6;
    double w = 2.0;
    double f_min = 0.01;
    double f_max = 2.0;
    int Nf = 512;
    int Ntheta = 32;
    int disable_math = 0;

    // Internal cache state
    bool initialized = false;
    bool dirty = true;

    std::vector<double> amps;
    std::vector<double> kx;
    std::vector<double> ky;
    std::vector<double> phis;

    // Cached parameters
    int last_seed = -1;
    int last_road_class = -1;
    double last_Gd_n0 = -1.0;
    double last_w = -1.0;
    double last_f_min = -1.0;
    double last_f_max = -1.0;
    int last_Nf = -1;
    int last_Ntheta = -1;

    void log(fmi2Status status, const char* category, const char* message) {
        if (logger) {
            logger(componentEnvironment, instanceName.c_str(), status, category, "%s", message);
        }
    }
};

// Numerical integration of (1 + t^2)^(-alpha/2)
double get_I(double alpha) {
    double sum = 0.0;
    double t_min = -2000.0;
    double t_max = 2000.0;
    int steps = 200000;
    double dt = (t_max - t_min) / steps;
    for (int i = 0; i < steps; ++i) {
        double t = t_min + i * dt;
        sum += std::pow(1.0 + t*t, -alpha/2.0);
    }
    return sum * dt;
}

void lazy_init(ModelInstance* inst) {
    if (inst->initialized &&
        !inst->dirty &&
        inst->last_seed == inst->seed &&
        inst->last_road_class == inst->road_class &&
        inst->last_Gd_n0 == inst->Gd_n0 &&
        inst->last_w == inst->w &&
        inst->last_f_min == inst->f_min &&
        inst->last_f_max == inst->f_max &&
        inst->last_Nf == inst->Nf &&
        inst->last_Ntheta == inst->Ntheta) {
        return;
    }

    NumPyRandom rng(inst->seed);

    // Map road class to Gd(n0)
    double Gd_n0_val = inst->Gd_n0;
    if (inst->road_class == 1) Gd_n0_val = 16e-6;
    else if (inst->road_class == 2) Gd_n0_val = 64e-6;
    else if (inst->road_class == 3) Gd_n0_val = 256e-6;
    else if (inst->road_class == 4) Gd_n0_val = 1024e-6;
    else if (inst->road_class == 5) Gd_n0_val = 4096e-6;

    double n0 = 0.1;
    double C1 = Gd_n0_val * std::pow(n0, inst->w);

    double alpha = inst->w + 1.0;
    double I_val = get_I(alpha);
    double C2 = C1 / I_val;

    int Nf = inst->Nf;
    int Ntheta = inst->Ntheta;

    // Logarithmic frequency spacing
    std::vector<double> f_r(Nf + 1);
    double log_f_min = std::log10(inst->f_min);
    double log_f_max = std::log10(inst->f_max);
    double dlog_f = (log_f_max - log_f_min) / Nf;
    for (int i = 0; i <= Nf; ++i) {
        f_r[i] = std::pow(10.0, log_f_min + i * dlog_f);
    }

    std::vector<double> df_r(Nf);
    std::vector<double> f_centers(Nf);
    for (int i = 0; i < Nf; ++i) {
        df_r[i] = f_r[i + 1] - f_r[i];
        f_centers[i] = 0.5 * (f_r[i] + f_r[i + 1]);
    }

    // Angular grid over [0, pi)
    double dtheta = M_PI / Ntheta;

    inst->amps.clear();
    inst->kx.clear();
    inst->ky.clear();
    inst->phis.clear();

    inst->amps.reserve(Nf * Ntheta);
    inst->kx.reserve(Nf * Ntheta);
    inst->ky.reserve(Nf * Ntheta);
    inst->phis.reserve(Nf * Ntheta);

    const double double_pi = 2.0 * M_PI;

    for (int i = 0; i < Nf; ++i) {
        double fc = f_centers[i];
        double dfc = df_r[i];
        double S_2D_val = C2 * std::pow(fc, -alpha);
        double power_per_angle = S_2D_val * fc * dfc * dtheta;
        double amp = std::sqrt(2.0 * power_per_angle);

        for (int j = 0; j < Ntheta; ++j) {
            double th = j * dtheta;
            double phi = rng.uniform(0.0, double_pi);

            inst->amps.push_back(amp);
            inst->kx.push_back(double_pi * fc * std::cos(th));
            inst->ky.push_back(double_pi * fc * std::sin(th));
            inst->phis.push_back(phi);
        }
    }

    // Cache the parameters
    inst->last_seed = inst->seed;
    inst->last_road_class = inst->road_class;
    inst->last_Gd_n0 = inst->Gd_n0;
    inst->last_w = inst->w;
    inst->last_f_min = inst->f_min;
    inst->last_f_max = inst->f_max;
    inst->last_Nf = inst->Nf;
    inst->last_Ntheta = inst->Ntheta;

    inst->initialized = true;
    inst->dirty = false;
}

inline double fast_cos(double x) {
    const double inv_two_pi = 0.15915494309189535;
    const double two_pi = 6.283185307179586;
    double z = x * inv_two_pi;
    double ip = std::floor(z + 0.5);
    double y = x - ip * two_pi;
    double y2 = y * y;
    return 1.0 + y2 * (-0.49999999682 + y2 * (0.041666641409 + y2 * (-0.0013888567281 + y2 * (0.000024786700184 + y2 * (-0.00000027246305604 + y2 * 0.0000000017860045536)))));
}

double compute_height(ModelInstance* inst, double px, double py) {
    if (inst->disable_math == 1) {
        return 0.0;
    }
    lazy_init(inst);

    double sum = 0.0;
    
    size_t size = inst->amps.size();
    const double* amps = inst->amps.data();
    const double* kx = inst->kx.data();
    const double* ky = inst->ky.data();
    const double* phis = inst->phis.data();

    for (int i = 0; i < (int)size; ++i) {
        double phase = kx[i] * px + ky[i] * py + phis[i];
        sum += amps[i] * fast_cos(phase);
    }
    return sum;
}

extern "C" {

FMI2_Export const char* fmi2GetTypesPlatform() {
    return fmi2TypesPlatform;
}

FMI2_Export const char* fmi2GetVersion() {
    return fmi2Version;
}

FMI2_Export fmi2Status fmi2SetDebugLogging(fmi2Component c, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    inst->loggingOn = loggingOn;
    return fmi2OK;
}

FMI2_Export fmi2Component fmi2Instantiate(fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGUID, fmi2String fmuResourceLocation, const fmi2CallbackFunctions* functions, fmi2Boolean visible, fmi2Boolean loggingOn) {
    if (fmuType != fmi2CoSimulation) {
        return nullptr;
    }
    
    ModelInstance* inst = new ModelInstance();
    inst->instanceName = instanceName ? instanceName : "";
    inst->guid = fmuGUID ? fmuGUID : "";
    if (functions) {
        inst->logger = functions->logger;
        inst->componentEnvironment = functions->componentEnvironment;
    }
    inst->loggingOn = loggingOn;
    inst->initialized = false;
    inst->dirty = true;
    
    return (fmi2Component)inst;
}

FMI2_Export void fmi2FreeInstance(fmi2Component c) {
    ModelInstance* inst = (ModelInstance*)c;
    if (inst) {
        delete inst;
    }
}

FMI2_Export fmi2Status fmi2SetupExperiment(fmi2Component c, fmi2Boolean toleranceDefined, fmi2Real tolerance, fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2EnterInitializationMode(fmi2Component c) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2ExitInitializationMode(fmi2Component c) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    lazy_init(inst);
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2Terminate(fmi2Component c) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2Reset(fmi2Component c) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    inst->initialized = false;
    inst->dirty = true;
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[]) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    for (size_t i = 0; i < nvr; ++i) {
        fmi2ValueReference v = vr[i];
        if (v == 0) {
            value[i] = inst->x;
        } else if (v == 1) {
            value[i] = inst->y;
        } else if (v == 2) {
            value[i] = compute_height(inst, inst->x, inst->y);
        } else if (v == 5) {
            value[i] = inst->Gd_n0;
        } else if (v == 6) {
            value[i] = inst->w;
        } else if (v == 7) {
            value[i] = inst->f_min;
        } else if (v == 8) {
            value[i] = inst->f_max;
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[]) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    for (size_t i = 0; i < nvr; ++i) {
        fmi2ValueReference v = vr[i];
        if (v == 0) {
            inst->x = value[i];
        } else if (v == 1) {
            inst->y = value[i];
        } else if (v == 2) {
            // z is read-only output, ignore writing
        } else if (v == 5) {
            inst->Gd_n0 = value[i];
            inst->dirty = true;
        } else if (v == 6) {
            inst->w = value[i];
            inst->dirty = true;
        } else if (v == 7) {
            inst->f_min = value[i];
            inst->dirty = true;
        } else if (v == 8) {
            inst->f_max = value[i];
            inst->dirty = true;
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[]) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    for (size_t i = 0; i < nvr; ++i) {
        fmi2ValueReference v = vr[i];
        if (v == 3) {
            value[i] = inst->seed;
        } else if (v == 4) {
            value[i] = inst->road_class;
        } else if (v == 9) {
            value[i] = inst->Nf;
        } else if (v == 10) {
            value[i] = inst->Ntheta;
        } else if (v == 11) {
            value[i] = inst->disable_math;
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[]) {
    ModelInstance* inst = (ModelInstance*)c;
    if (!inst) return fmi2Error;
    for (size_t i = 0; i < nvr; ++i) {
        fmi2ValueReference v = vr[i];
        if (v == 3) {
            inst->seed = value[i];
            inst->dirty = true;
        } else if (v == 4) {
            inst->road_class = value[i];
            inst->dirty = true;
        } else if (v == 9) {
            inst->Nf = value[i];
            inst->dirty = true;
        } else if (v == 10) {
            inst->Ntheta = value[i];
            inst->dirty = true;
        } else if (v == 11) {
            inst->disable_math = value[i];
        } else {
            return fmi2Error;
        }
    }
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[]) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[]) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[]) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[]) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2DoStep(fmi2Component c, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize, fmi2Boolean noSetFMUStatePriorToCurrentPoint) {
    return fmi2OK;
}

FMI2_Export fmi2Status fmi2CancelStep(fmi2Component c) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind s, fmi2Status* value) {
    return fmi2Discard;
}

FMI2_Export fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind s, fmi2Real* value) {
    return fmi2Discard;
}

FMI2_Export fmi2Status fmi2GetIntegerStatus(fmi2Component c, const fmi2StatusKind s, fmi2Integer* value) {
    return fmi2Discard;
}

FMI2_Export fmi2Status fmi2GetBooleanStatus(fmi2Component c, const fmi2StatusKind s, fmi2Boolean* value) {
    return fmi2Discard;
}

FMI2_Export fmi2Status fmi2GetStringStatus(fmi2Component c, const fmi2StatusKind s, fmi2String* value) {
    return fmi2Discard;
}

FMI2_Export fmi2Status fmi2GetFMUstate(fmi2Component c, fmi2FMUstate* FMUstate) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2SetFMUstate(fmi2Component c, fmi2FMUstate FMUstate) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* FMUstate) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate FMUstate, size_t* size) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2SerializeFMUstate(fmi2Component c, fmi2FMUstate FMUstate, fmi2Byte serializedState[], size_t size) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2DeSerializeFMUstate(fmi2Component c, const fmi2Byte serializedState[], size_t size, fmi2FMUstate* FMUstate) {
    return fmi2Error;
}

FMI2_Export fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vUnknown[], size_t nUnknown, const fmi2ValueReference vKnown[], size_t nKnown, const fmi2Real dvKnown[], fmi2Real dvUnknown[]) {
    return fmi2Error;
}

// Model Exchange Dummy Functions to satisfy full linking requirements
FMI2_Export fmi2Status fmi2EnterEventMode(fmi2Component c) { return fmi2Error; }
FMI2_Export fmi2Status fmi2NewDiscreteStates(fmi2Component c, fmi2EventInfo* eventInfo) { return fmi2Error; }
FMI2_Export fmi2Status fmi2EnterContinuousTimeMode(fmi2Component c) { return fmi2Error; }
FMI2_Export fmi2Status fmi2CompletedIntegratorStep(fmi2Component c, fmi2Boolean noSetFMUStatePriorToCurrentPoint, fmi2Boolean* enterEventMode, fmi2Boolean* terminateSimulation) { return fmi2Error; }
FMI2_Export fmi2Status fmi2SetTime(fmi2Component c, fmi2Real time) { return fmi2Error; }
FMI2_Export fmi2Status fmi2SetContinuousStates(fmi2Component c, const fmi2Real x[], size_t nx) { return fmi2Error; }
FMI2_Export fmi2Status fmi2GetDerivatives(fmi2Component c, fmi2Real derivatives[], size_t nx) { return fmi2Error; }
FMI2_Export fmi2Status fmi2GetEventIndicators(fmi2Component c, fmi2Real eventIndicators[], size_t ni) { return fmi2Error; }
FMI2_Export fmi2Status fmi2GetContinuousStates(fmi2Component c, fmi2Real x[], size_t nx) { return fmi2Error; }
FMI2_Export fmi2Status fmi2GetNominalsOfContinuousStates(fmi2Component c, fmi2Real x_nominal[], size_t nx) { return fmi2Error; }
FMI2_Export fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], const fmi2Real value[]) { return fmi2Error; }
FMI2_Export fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], fmi2Real value[]) { return fmi2Error; }

} // extern "C"
