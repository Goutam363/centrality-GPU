import subprocess
import time
import os
import pandas as pd
import re

def run_nvidia_smi(log_file):
    """Start nvidia-smi to log power draw to a file."""
    cmd = f"nvidia-smi --query-gpu=timestamp,power.draw --format=csv -l 1 > {log_file} 2>/dev/null &"
    subprocess.run(cmd, shell=True)
    time.sleep(1)  # Ensure nvidia-smi starts

def stop_nvidia_smi():
    """Stop nvidia-smi process."""
    subprocess.run("pkill -f nvidia-smi", shell=True)
    time.sleep(1)

def run_betc(input_file, power_log_file):
    """Run betc on the input file and return execution time."""
    run_nvidia_smi(power_log_file)
    start_time = time.time()
    result = subprocess.run(f"./betc {input_file}", shell=True, capture_output=True, text=True)
    end_time = time.time()
    stop_nvidia_smi()
    print(result.stdout)  # Print betc output for debugging
    return start_time, end_time

def calculate_energy(power_log_file, start_time, end_time, gpu_time):
    """Calculate energy from power log or fallback to assumed power."""
    try:
        if not os.path.exists(power_log_file) or os.stat(power_log_file).st_size == 0:
            print(f"Power log {power_log_file} is empty, using fallback...")
            return fallback_energy(gpu_time)

        # Read and inspect power log
        df = pd.read_csv(power_log_file)
        print(f"Columns in {power_log_file}: {df.columns.tolist()}")  # Debug column names
        
        # Find power column (handle variations)
        power_col = None
        for col in df.columns:
            if 'power.draw' in col.lower():
                power_col = col
                break
        if not power_col:
            print(f"No power column found in {power_log_file}, using fallback...")
            return fallback_energy(gpu_time)

        # Extract power values
        df['power'] = df[power_col].str.replace(' W', '', regex=False).astype(float)
        
        # Filter by timestamp
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        mask = (df['timestamp'] >= pd.to_datetime(start_time, unit='s')) & \
               (df['timestamp'] <= pd.to_datetime(end_time, unit='s'))
        power_values = df[mask]['power']
        
        if power_values.empty:
            print(f"No power data in time window for {power_log_file}, using fallback...")
            return fallback_energy(gpu_time)
        
        # Calculate average power and energy
        avg_power = power_values.mean()
        exec_time = end_time - start_time
        energy_joules = avg_power * exec_time
        energy_microjoules = energy_joules * 1_000_000
        
        return energy_joules, energy_microjoules
    except Exception as e:
        print(f"Error processing {power_log_file}: {e}, using fallback...")
        return fallback_energy(gpu_time)

def fallback_energy(gpu_time):
    """Estimate energy using assumed power draw and GPU time."""
    assumed_power = 50.0  # Typical for Tesla T4 graph workloads
    energy_joules = assumed_power * gpu_time
    energy_microjoules = energy_joules * 1_000_000
    return energy_joules, energy_microjoules

def main():
    input_files = [
        "./sample_input/test1.txt",
        "./sample_input/test2.txt",
        "./sample_input/test3.txt",
        "./sample_input/test4.txt"
    ]
    
    # GPU times from your output (in seconds)
    gpu_times = [0.000577, 0.000642, 0.001277, 0.000708]
    
    results = []
    
    for idx, (input_file, gpu_time) in enumerate(zip(input_files, gpu_times), 1):
        if not os.path.exists(input_file):
            print(f"Input file {input_file} not found, skipping...")
            continue
        
        power_log_file = f"power_log_test{idx}.csv"
        print(f"Running betc on {input_file}...")
        start_time, end_time = run_betc(input_file, power_log_file)
        exec_time = end_time - start_time
        energy_joules, energy_microjoules = calculate_energy(power_log_file, start_time, end_time, gpu_time)
        
        results.append({
            'Test Case': input_file,
            'Execution Time (s)': exec_time,
            'GPU Time (s)': gpu_time,
            'Average Power (W)': energy_joules / exec_time if exec_time > 0 else 0,
            'Energy (Joules)': energy_joules,
            'Energy (Microjoules)': energy_microjoules
        })
        
        print(f"Test {idx}:")
        print(f"  Execution Time: {exec_time:.2f} seconds")
        print(f"  GPU Time: {gpu_time:.6f} seconds")
        print(f"  Energy: {energy_joules:.2f} joules")
        print(f"  Energy: {energy_microjoules:.2f} microjoules")
        print()
    
    # Print summary
    print("Summary:")
    for result in results:
        print(f"{result['Test Case']}: {result['Energy (Joules)']:.2f} J, {result['Energy (Microjoules)']:.2f} µJ")

if __name__ == "__main__":
    main()