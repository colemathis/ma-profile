import pandas as pd
import molecular_assembly as ma

INPUTNAME = "KEGG_complete_conversion_26052022.csv"
OUTPUTNAME = "KEGG_time_est.csv"
TIMEOUT = 60
def main():

    input_df = pd.read_csv(INPUTNAME)
    all_inchis = list(input_df["inchi"])
    time_est_maps = dict()
    n_inchis = len(all_inchis)
    print(n_inchis)
    for i in range(5):
        inchi = all_inchis[i]
        if i % 100 == 1:
            print(i)
        checked_inchi, covalent = ma.check_inchi(inchi, strict=False)
        size = ma.mol_and_inchis.check_small(checked_inchi)
        if size > 3:
            # print(checked_inchi)
            this_est = ma.estimate_time(checked_inchi, timeout=TIMEOUT)
            time_est_maps[inchi] = this_est
            # print(this_est)

    input_df["time_est"] = input_df["inchi"].map(time_est_maps)
    input_df.to_csv(OUTPUTNAME)


if __name__ == "__main__":
    main()
