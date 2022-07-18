import molecular_assembly as ma
import pandas as pd

SAMPLE_FNAME = "Sample_Benchmark_3-26.csv"
SAVENAME = "Sample_Subgraphs.csv"
TIMEOUT = 200

def main():
    """Compute the Number of Subgraphs in the benchmark compounds."""
    source_df = pd.read_csv("Sample_Benchmark_3-26.csv")
    inchis = list(source_df["InChI"])
    all_data = []

    for inchi in inchis:
        checked_inchi, _ = ma.check_inchi(inchi, strict=False)
        ma_run = ma.get_subgraphs(checked_inchi,
                                  timeout=TIMEOUT)
        save_data = {"subgraphs": ma_run["subgraphs"],
                     "InChI": inchi}
        all_data.append(save_data)

    time_df = pd.DataFrame(all_data)
    time_df.to_csv(SAVENAME)


if __name__ == "__main__":
    main()
