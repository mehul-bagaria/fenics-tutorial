import subprocess
import os

def run_fenics_simulation():
    """Runs the FEniCS simulation inside a Docker container."""
    current_dir = os.getcwd()
    print(f"Current directory: {current_dir}")
    
    docker_command = [
        "docker", "run", 
        "--rm",
        "-v", f"{current_dir}:/root",
        "-w", "/root",
        "iitrabhi/fenics",
        "python3 run_simulation.py"
    ]
    
    print("Starting FEniCS simulation...")
    
    try:
        result = subprocess.run(docker_command)
        
        if result.returncode == 0:
            print(" Simulation completed successfully!")
        else:
            print(f" Simulation failed with exit code: {result.returncode}")
        
        return result.returncode
        
    except Exception as e:
        print(f" Error: {e}")
        return 1

if __name__ == "__main__":
    exit(run_fenics_simulation())
