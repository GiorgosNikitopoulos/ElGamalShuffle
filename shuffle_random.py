import random
import switch
from random import randint

RANDOM_FUNC_CHOICE = 1
def shuffle_rand_int(options, floor, ceiling, **addition_arguments):
    """Wrapper for generating pseudorandom numbers"""
    """Index: 1: Using the randint function from the random module
              2: Using the more cryptographically secure secrets module and
              the rand below function
    """
    while switch.switch(options):
        if switch.case(1):
            return randint(floor, ceiling)
            break
        if switch.case(2):
            rng = random.SystemRandom()
            return rng.random(floor, ceiling)
            break

        # Add cases similarly and call functions appropriately.
        # Options changes the pseudorandom generation method.
        # Use addition_arguments for any additional arguments needed for the operation
        # of any pseudorandom generation function that needs more arguments.
        break
