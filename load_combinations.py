# Load cases (only DL and LL as requested)
load_cases = ["DL", "LL"]

# Load combinations
load_combinations = {
    "STR_1": {"DL": 1.4, "LL": 1.6},
    "SER_1": {"DL": 1.0, "LL": 1.0},
    "STR_2": {"DL": 1.2, "LL": 1.4}
}

# Node names
node_names = ["N1", "N2", "N3"]

# Manually created json_data
json_data = [
    # Node N1
    {
        
        "node_name": "N1",
        "load_combination": "STR_1",
        "load_case": "DL",
        "load_values": [120.5, 60.2, 1500.0, 15.3, 30.1, 7.5],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N1",
        "load_combination": "STR_1",
        "load_case": "LL",
        "load_values": [80.3, 40.1, 1000.0, 10.2, 20.0, 5.0],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N1",
        "load_combination": "SER_1",
        "load_case": "DL",
        "load_values": [86.1, 43.0, 1071.4, 10.9, 21.5, 5.4],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N1",
        "load_combination": "SER_1",
        "load_case": "LL",
        "load_values": [80.3, 40.1, 1000.0, 10.2, 20.0, 5.0],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    
    # Node N2
    {
        "node_name": "N2",
        "load_combination": "STR_1",
        "load_case": "DL",
        "load_values": [110.2, 55.1, 1375.0, 13.8, 27.5, 6.9],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N2",
        "load_combination": "STR_1",
        "load_case": "LL",
        "load_values": [75.0, 37.5, 937.5, 9.4, 18.8, 4.7],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N2",
        "load_combination": "STR_2",
        "load_case": "DL",
        "load_values": [94.4, 47.2, 1180.0, 11.8, 23.6, 5.9],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N2",
        "load_combination": "STR_2",
        "load_case": "LL",
        "load_values": [105.0, 52.5, 1312.5, 13.1, 26.3, 6.6],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    
    # Node N3
    {
        "node_name": "N3",
        "load_combination": "SER_1",
        "load_case": "DL",
        "load_values": [95.0, 47.5, 1187.5, 11.9, 23.8, 5.9],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N3",
        "load_combination": "SER_1",
        "load_case": "LL",
        "load_values": [85.0, 42.5, 1062.5, 10.6, 21.3, 5.3],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N3",
        "load_combination": "STR_2",
        "load_case": "DL",
        "load_values": [114.0, 57.0, 1425.0, 14.3, 28.5, 7.1],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    },
    {
        "node_name": "N3",
        "load_combination": "STR_2",
        "load_case": "LL",
        "load_values": [119.0, 59.5, 1487.5, 14.9, 29.8, 7.4],
        "units": {
            "FX": "kN", "FY": "kN", "FZ": "kN",
            "MX": "kN-m", "MY": "kN-m", "MZ": "kN-m"
        }
    }
]
