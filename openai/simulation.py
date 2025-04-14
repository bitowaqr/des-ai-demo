import random
import math
import json
import numpy as np

def simulate_patient(treatment, params):
    """
    Simulates one patient’s trajectory.
    
    For each patient, we sample:
      - t_stroke: time to stroke event (based on treatment-specific hazard)
      - t_death: time to death (using a natural annual mortality rate)
    
    If the stroke happens before death (and before the simulation horizon),
    we split the timeline into pre-stroke and post-stroke periods with different 
    utility and cost implications. Continuous discounting is applied.
    """
    horizon = params["horizon"]
    r = params["discount_rate"]
    
    # Choose treatment-specific parameters
    if treatment == "supimab":
        stroke_rate = params["baseline_stroke_rate"] * params["hazard_ratio_supimab"]
        drug_cost = params["drug_cost_supimab"]
    else:
        stroke_rate = params["baseline_stroke_rate"]
        drug_cost = params["drug_cost_dupimab"]
    
    t_stroke = random.expovariate(stroke_rate)
    t_death = random.expovariate(params["mortality_rate_natural"])
    
    # Determine if stroke occurs (only if it happens before death and within horizon)
    if t_death < t_stroke or t_stroke > horizon:
        t_end = min(t_death, horizon)
        stroke_occurred = False
        # Calculate QALYs and treatment cost over the whole period using the pre-stroke utility
        qalys = params["utility_pre_stroke"] / r * (1 - math.exp(-r * t_end))
        cost = drug_cost / r * (1 - math.exp(-r * t_end))
    else:
        # Stroke event occurs before death and within horizon
        t_end = min(t_death, horizon)
        stroke_occurred = True
        # Pre-stroke period (0 to t_stroke)
        qalys_pre = params["utility_pre_stroke"] / r * (1 - math.exp(-r * t_stroke))
        cost_pre = drug_cost / r * (1 - math.exp(-r * t_stroke))
        # Immediate stroke event cost (discounted to time t_stroke)
        stroke_cost_disc = params["stroke_cost"] * math.exp(-r * t_stroke)
        # Post-stroke period (from t_stroke to t_end)
        qalys_post = params["utility_post_stroke"] / r * (math.exp(-r * t_stroke) - math.exp(-r * t_end))
        cost_post_drug = drug_cost / r * (math.exp(-r * t_stroke) - math.exp(-r * t_end))
        cost_post_stroke = params["annual_poststroke_cost"] / r * (math.exp(-r * t_stroke) - math.exp(-r * t_end))
        qalys = qalys_pre + qalys_post
        cost = cost_pre + stroke_cost_disc + cost_post_drug + cost_post_stroke

    return {
        "treatment": treatment,
        "cost": cost,
        "qalys": qalys,
        "stroke_occurred": stroke_occurred,
        "t_stroke": t_stroke if stroke_occurred else None,
        "t_death": t_death,
        "t_end": t_end
    }

def run_simulation(params):
    """Simulate a cohort for each treatment arm."""
    results = {"dupimab": [], "supimab": []}
    for treatment in ["dupimab", "supimab"]:
        for i in range(params["n_patients"]):
            results[treatment].append(simulate_patient(treatment, params))
    return results

def summarize_results(sim_data, treatment):
    """
    Computes summary statistics (mean, standard deviation) for cost and QALYs.
    Also returns individual patient outcomes for further analyses.
    """
    costs = [p["cost"] for p in sim_data[treatment]]
    qalys = [p["qalys"] for p in sim_data[treatment]]
    mean_cost = np.mean(costs)
    mean_qalys = np.mean(qalys)
    std_cost = np.std(costs)
    std_qalys = np.std(qalys)
    return {
        "mean_cost": mean_cost,
        "mean_qalys": mean_qalys,
        "std_cost": std_cost,
        "std_qalys": std_qalys,
        "costs": costs,
        "qalys": qalys
    }

def main():
    # Model parameters (modifiable)
    params = {
        "horizon": 30,                    # Simulation horizon in years
        "discount_rate": 0.035,           # Continuous discount rate (3.5%)
        "n_patients": 1000,               # Number of simulated patients per arm
        "baseline_stroke_rate": 0.05,     # Annual stroke risk for dupimab
        "hazard_ratio_supimab": 0.7,      # Supimab reduces stroke risk by 30%
        "utility_pre_stroke": 0.85,       # Utility before stroke
        "utility_post_stroke": 0.70,      # Utility after stroke
        "drug_cost_dupimab": 10000,       # Annual cost for dupimab (£)
        "drug_cost_supimab": 15000,       # Annual cost for supimab (£)
        "stroke_cost": 20000,             # One-off cost at stroke event (£)
        "annual_poststroke_cost": 5000,   # Additional annual cost after stroke (£)
        "mortality_rate_natural": 0.02    # Annual natural mortality rate
    }
    
    sim_data = run_simulation(params)
    
    summary_dupimab = summarize_results(sim_data, "dupimab")
    summary_supimab = summarize_results(sim_data, "supimab")
    
    # Compute the incremental cost-effectiveness ratio (ICER)
    icer = (summary_supimab["mean_cost"] - summary_dupimab["mean_cost"]) / (summary_supimab["mean_qalys"] - summary_dupimab["mean_qalys"])
    
    output = {
        "parameters": params,
        "summary": {
            "dupimab": summary_dupimab,
            "supimab": summary_supimab,
            "icer": icer
        },
        "individuals": sim_data
    }
    
    # Write simulation results to JSON for the dashboard
    with open("simulation_results.json", "w") as f:
        json.dump(output, f, indent=4)
    print("Simulation complete. Results saved to simulation_results.json")

if __name__ == "__main__":
    main()

