import sys
import importlib
import subprocess
import shutil
import os

def get_tool_version(command, version_flag='--version'):
    try:
        executable = command.split()[0]
        if executable != 'java' and shutil.which(executable) is None and not os.path.exists(executable):
             return "not_found"
        
        try:
            output = subprocess.check_output(f"{command} {version_flag}", shell=True, stderr=subprocess.STDOUT).decode().strip()
        except subprocess.CalledProcessError as e:
            output = e.output.decode().strip()
            
        lines = output.split('\n')
        
        if 'bwa-mem2' in command:
            for line in lines:
                if "Looking to launch" in line:
                    continue
                if "Version:" in line:
                    return line.split('Version:')[1].strip()
                if any(c.isdigit() for c in line) and '.' in line:
                    return line.strip()

        if 'gatk' in command:
             for line in lines:
                if 'The GATK version is' in line:
                    return line.split('The GATK version is')[1].strip()

        for line in lines:
            if "Looking to launch" in line:
                continue
            if any(c.isdigit() for c in line):
                return line.strip().replace('"', '')
                
        return lines[0].strip()
    except Exception as e:
        return f"error"

def main():
    packages = [
        "pandas", "numpy"
    ]

    print("python_tools:")
    print(f"  python: {sys.version.split()[0]}")
    
    for p in packages:
        try:
            mod = importlib.import_module(p)
            version = getattr(mod, "__version__", "unknown")
            print(f"  {p}: {version}")
        except ImportError:
            print(f"  {p}: not_installed")

    print("system_tools:")
    
    base_dir = os.environ.get('BASE_DIR', '.')
    gatk_path = os.path.join(base_dir, "tools/gatk-4.4.0.0/gatk")
    
    tools = {
        "bcftools": ("bcftools", "--version"),
        "samtools": ("samtools", "--version"),
        "bwa-mem2": ("bwa-mem2", "version"),
        "minimap2": ("minimap2", "--version"),
        "fastqc": ("fastqc", "--version"),
        "multiqc": ("multiqc", "--version"),
        "java": ("java", "-version"),
        "medaka": ("medaka", "--version"),
        "chopper": ("chopper", "--version"),
        "gatk": (gatk_path, "--version"),
        "trimmomatic": ("java -jar /usr/share/java/trimmomatic.jar", "-version")
    }

    for name, (cmd, flag) in tools.items():
        ver = get_tool_version(cmd, flag)
        print(f"  {name}: \"{ver}\"")

if __name__ == "__main__":
    main()