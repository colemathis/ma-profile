import molecular_assembly as ma
import pandas as pd

TIMEOUT = 1000
SAVENAME = "Exact_Go_Times_1000s_2022_07_14.csv"


def main():
    """Compute the MA of the benchmark compounds and record the time."""
    source_df = pd.read_csv("Sample_Benchmark_3-26.csv")
    inchis = list(source_df["InChI"])
    all_data = []
    for inchi in inchis:
        # checked_inchi, _ = ma.check_inchi(inchi, strict=False)
        # print(checked_inchi)
        ma_run = ma.get_exact_go(inchi, timeout=TIMEOUT)
        saved_data = {"MA": ma_run["MA"],
                      "InChI": inchi,
                      "time": ma_run["time"]}
        all_data.append(saved_data)

    time_df = pd.DataFrame(all_data)
    time_df.to_csv(SAVENAME)


if __name__ == "__main__":
    main()
