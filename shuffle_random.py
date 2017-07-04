import random
import switch
from random import randint


def shuffle_rand_int(options, floor, ceiling, **addition_arguments):
    """Wrapper for generating pseudorandom numbers"""
    while switch.switch(options):
        if switch.case(1):
            return randint(floor, ceiling)
            break
        # Add cases similarly and call functions appropriately.
        # Options changes the pseudorandom generation method.
        # Use addition_arguments for any additional arguments needed for the operation
        # of any pseudorandom generation function that needs more arguments.
        break
