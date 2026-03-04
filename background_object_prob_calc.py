""" 
Background object probability calculator
MPhys Project 
Aiza Kenzhebekova 
02/03/2026
"""
import background_object_prob_calc_functions as bckgprob 
import argparse

def main():
    """ 
    This calculates the probability of finding a background object w/ specified parameters
    ----------------------------------------------------------------------------------------
    Format: python background_object_calc.py "inputfile" --trilegal RA Dec FOV
    TRILEGAL run is optional. If the colours of your objects are red, skip trilegal. 

    Inputs: 
    "inputfile" = .dat file where the header is either "F444W F200W-F444W Separation" 
            or "F1500W F1500W-F2100W Separation" 
    RA, Dec in degrees
    FOV in sq. arcsec

    Output: Probability (%)
    """
    parser = argparse.ArgumentParser(description="Calculate background object probability.")
    
    # required argument (input file)
    parser.add_argument("input_file", help="Input data file (.dat) ")

    # optional arguments (ra dec fov)
    parser.add_argument("--trilegal", nargs=3, metavar=("RA", "DEC", "FOV"), type=float, 
                help="Run TRILEGAL for specific RA (deg), DEC (deg), and with a FOV (sq. arcsec)")
    args = parser.parse_args()

    # w/trilegal
    if args.trilegal:
        # 
        ra, dec, fov = args.trilegal
        probability = bckgprob.calculate_probability(args.input_file, 
                                 run_trilegal=True, trilegal_los=[ra, dec], trilegal_fov=fov)

    # w/o trilegal
    else:
        probability = bckgprob.calculate_probability(args.input_file)

    for i, prob in enumerate(probability):
        print(f"Row {i+1}: {round(prob, 4)}%")

if __name__ == "__main__":
    main()