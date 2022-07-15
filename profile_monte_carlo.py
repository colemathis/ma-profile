import molecular_assembly as ma
import pandas as pd

SAVENAME = "Monte_Carlo_Times_200s_2022_07_14.csv"
TIMEOUT = 10
MIN_MA = 14
FIRST_STEPS = [1000, 10000, 100000]
SECOND_STEPS = [100, 500, 1000, int(1e6)]


def main():
    """Compute the MA of the benchmark compounds and record the time."""
    source_df = pd.read_csv("Sample_Benchmark_3-26.csv")
    source_df = source_df[source_df["MA"] > MIN_MA]
    inchis = list(source_df["InChI"])
    all_data = []

    for nstep1 in FIRST_STEPS:
        for nstep2 in SECOND_STEPS:
            for inchi in inchis[:2]:
                checked_inchi = ma.check_inchi_c(inchi, strict=False)
                ma_run = ma.get_monte_carlo(checked_inchi,
                                            Nstep1=nstep1,
                                            Nstep2=nstep2,
                                            timeout=TIMEOUT)
                save_data = {"MA": ma_run["MA"],
                             "time": ma_run["time"],
                             "InChI": inchi,
                             "Nstep1": nstep1,
                             "Nstep2": nstep2}
                all_data.append(save_data)

    time_df = pd.DataFrame(all_data)
    time_df.to_csv(SAVENAME)


if __name__ == "__main__":
    main()
