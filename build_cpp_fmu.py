import os
import sys
import shutil
import subprocess
import zipfile

def build_fmu():
    print("=== STARTING C++ FMU BUILD & PACKAGING ===")
    
    # 1. Compile C++ Shared Library using CMake
    src_dir = os.path.abspath("cpp_fmu")
    build_dir = os.path.join(src_dir, "build")
    
    print(f"Configuring CMake in {build_dir}...")
    os.makedirs(build_dir, exist_ok=True)
    configure_cmd = ["cmake", "-B", "build", "-S", "."]
    res = subprocess.run(configure_cmd, cwd=src_dir)
    if res.returncode != 0:
        print("Error: CMake configuration failed.")
        sys.exit(1)
        
    print("Compiling C++ shared library in Release mode...")
    build_cmd = ["cmake", "--build", "build", "--config", "Release"]
    res = subprocess.run(build_cmd, cwd=src_dir)
    if res.returncode != 0:
        print("Error: Compilation failed.")
        sys.exit(1)
        
    # Check compiled DLL
    dll_source = os.path.join(build_dir, "Release", "InfiniteRoadFMU.dll")
    if not os.path.exists(dll_source):
        # Check direct build dir just in case
        dll_source = os.path.join(build_dir, "InfiniteRoadFMU.dll")
        
    if not os.path.exists(dll_source):
        print(f"Error: Compiled DLL not found at expected path: {dll_source}")
        sys.exit(1)
        
    print(f"Compiled DLL found at: {dll_source}")
    
    # 2. Setup Staging Directory for FMU packaging
    staging_dir = os.path.abspath("fmu_staging")
    if os.path.exists(staging_dir):
        shutil.rmtree(staging_dir)
    os.makedirs(staging_dir)
    
    bin_win64_dir = os.path.join(staging_dir, "binaries", "win64")
    os.makedirs(bin_win64_dir)
    
    # Copy DLL to FMI binaries folder
    dll_dest = os.path.join(bin_win64_dir, "InfiniteRoadFMU.dll")
    shutil.copy2(dll_source, dll_dest)
    print(f"Staged DLL to {dll_dest}")
    
    # 3. Create modelDescription.xml in staging
    xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<fmiModelDescription fmiVersion="2.0" modelName="InfiniteRoadFMU" guid="ab782f63-6d6a-11f1-a954-8ce9ee8a1429" generationTool="Native C++ Implementation" generationDateAndTime="2026-06-21T15:42:00Z" variableNamingConvention="structured" description="Deterministic ISO 8608 2D Isotropic Infinite Road Profile Generator (Single Point)" author="Antigravity Coding Assistant">
	<CoSimulation needsExecutionTool="false" canHandleVariableCommunicationStepSize="true" canInterpolateInputs="false" canBeInstantiatedOnlyOncePerProcess="false" canGetAndSetFMUstate="false" canSerializeFMUstate="false" modelIdentifier="InfiniteRoadFMU" canNotUseMemoryManagementFunctions="true"/>
	<LogCategories>
		<Category name="logStatusWarning" description="Log messages with fmi2Warning status."/>
		<Category name="logStatusDiscard" description="Log messages with fmi2Discard status."/>
		<Category name="logStatusError" description="Log messages with fmi2Error status."/>
		<Category name="logStatusFatal" description="Log messages with fmi2Fatal status."/>
		<Category name="logAll" description="Log all messages."/>
	</LogCategories>
	<ModelVariables>
		<ScalarVariable name="x" valueReference="0" causality="input">
			<Real start="0"/>
		</ScalarVariable>
		<ScalarVariable name="y" valueReference="1" causality="input">
			<Real start="0"/>
		</ScalarVariable>
		<ScalarVariable name="z" valueReference="2" causality="output">
			<Real/>
		</ScalarVariable>
		<ScalarVariable name="seed" valueReference="3" causality="parameter" variability="fixed">
			<Integer start="42"/>
		</ScalarVariable>
		<ScalarVariable name="road_class" valueReference="4" causality="parameter" variability="fixed">
			<Integer start="3"/>
		</ScalarVariable>
		<ScalarVariable name="Gd_n0" valueReference="5" causality="parameter" variability="fixed">
			<Real start="0.000256"/>
		</ScalarVariable>
		<ScalarVariable name="w" valueReference="6" causality="parameter" variability="fixed">
			<Real start="2"/>
		</ScalarVariable>
		<ScalarVariable name="f_min" valueReference="7" causality="parameter" variability="fixed">
			<Real start="0.002"/>
		</ScalarVariable>
		<ScalarVariable name="f_max" valueReference="8" causality="parameter" variability="fixed">
			<Real start="2000"/>
		</ScalarVariable>
		<ScalarVariable name="Nf" valueReference="9" causality="parameter" variability="fixed">
			<Integer start="512"/>
		</ScalarVariable>
		<ScalarVariable name="Ntheta" valueReference="10" causality="parameter" variability="fixed">
			<Integer start="32"/>
		</ScalarVariable>
		<ScalarVariable name="disable_math" valueReference="11" causality="parameter" variability="fixed">
			<Integer start="0"/>
		</ScalarVariable>
	</ModelVariables>
	<ModelStructure>
		<Outputs>
			<Unknown index="3"/>
		</Outputs>
	</ModelStructure>
</fmiModelDescription>
"""
    xml_dest = os.path.join(staging_dir, "modelDescription.xml")
    with open(xml_dest, "w", encoding="utf-8") as f:
        f.write(xml_content.strip())
    print(f"Created modelDescription.xml at {xml_dest}")
    
    # 4. Zip the staging folder contents to create the final FMU
    fmu_filename = "InfiniteRoadFMU.fmu"
    if os.path.exists(fmu_filename):
        os.remove(fmu_filename)
        
    print(f"Creating zip archive {fmu_filename}...")
    with zipfile.ZipFile(fmu_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(staging_dir):
            for file in files:
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, staging_dir)
                zipf.write(abs_path, rel_path)
                
    # Clean up staging dir
    shutil.rmtree(staging_dir)
    print(f"Successfully packaged and created {fmu_filename}!")

if __name__ == "__main__":
    build_fmu()
