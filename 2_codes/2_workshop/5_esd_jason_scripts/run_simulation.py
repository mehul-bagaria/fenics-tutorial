import json
from src.static import run_plane_stress_simulation

if __name__ == "__main__":
    # Run the plane stress simulation with the configuration file
    config_path = "simulation_config/study_config.json"
    
    print(f"Running plane stress simulation with config: {config_path}")
    summary = run_plane_stress_simulation(config_path)
    
    # Print the summary results
    print("\n" + "="*60)
    print("SIMULATION SUMMARY")
    print("="*60)
    print(json.dumps(summary, indent=2))
    print("="*60)