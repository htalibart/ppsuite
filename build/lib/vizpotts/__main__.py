import argparse

from vizpotts import *
from potts_model import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--potts_model', help="Potts model (msgpack file)")
    #parser.add_argument('-m', '--mode', help="Visualization mode", choices=("whole_potts_model"), default="whole_potts_model")
    parser.add_argument('-1', '--start_at_1', help="Start numbering at 1", action='store_true', default=True), 
    parser.add_argument('-0', '--start_at_0', help="Start numbering at 0", action='store_true', default=False), 
    args = vars(parser.parse_args())


    start_at_1 = args["start_at_1"] and not args["start_at_0"]

    mrf = Potts_Model.from_msgpack(args["potts_model"])
    visualize_mrf(mrf, start_at_1=start_at_1)
