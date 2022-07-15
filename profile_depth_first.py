import molecular_assembly as ma
import pandas as pd

TIMEOUT = 10
SAVENAME = "Depth_First_Times_200s_2022_07_14.csv"


def main():
    """Compute the MA of the benchmark compounds and record the time."""
    source_df = pd.read_csv("Sample_Benchmark_3-26.csv")
    inchis = list(source_df["InChI"])
    all_data = []
    for inchi in inchis[50:55]:
        checked_inchi = ma.check_inchi_c(inchi, strict=False)
        ma_run = ma.get_depth_first(checked_inchi, timeout=TIMEOUT)
        save_data = {"MA": ma_run["MA"],
                     "time": ma_run["time"],
                     "InChI": inchi}
        all_data.append(save_data)

    time_df = pd.DataFrame(all_data)
    time_df.to_csv(SAVENAME)


if __name__ == "__main__":
    main()
