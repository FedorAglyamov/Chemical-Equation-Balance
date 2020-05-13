# Simple program for balancing chemical equations
# v0.9


# Imports
import numpy as np
from scipy import linalg
from fractions import Fraction
import re
import sys


# Declare constants
SEPERATOR = "="
TUTORIAL_MSG = "Please seperate left/right side of equation with '=' char, " \
               "and expand components like (OH)3 to O3H3."
BARS = "\n-------------------------------------------------------------------"
ATTEMPT_BAL_MSG = "\n*** Attempted Balance with {}***"


# Main function running program
def main():

    print(TUTORIAL_MSG)
    
    while True:
        print(BARS)
        # Find balanced equation
        reaction = get_reaction()
        coeff_matrix = build_coeff_matrix(reaction)
        print(coeff_matrix)
        frac_coeffs, coeffs = solve_coeff(coeff_matrix)
        print(ATTEMPT_BAL_MSG.format("Fraction Coefficients"))
        print_eq(reaction, frac_coeffs)
        print(ATTEMPT_BAL_MSG.format("Whole Coefficients"))
        print_eq(reaction, coeffs)
        # Confirm whether user would like to balance more equations
        if (not more_eq()):
            break
    

# Get user's input for chem equation and do some basic error checking
def get_reaction():

    # Initialize vars
    user_input = direction = ""
    reaction = []
    sides = []
    cur_side = []
    cur_elems = []
    side_elems = []
    elems = []

    # While user equation does not have two sides
    while True:
        user_input = input("\nPlease input reaction: ")
        # If user passes an equation
        if (user_input.count(SEPERATOR) == 1):
            # Get reaction and sides
            reaction = list(filter(None, re.split(r"\s|\+|(" + SEPERATOR + ")", user_input)))
            direction = re.findall(r"" + SEPERATOR, user_input)
            sides = list(filter(None, re.split(r"" + SEPERATOR, user_input)))
            # Find elems present on left and right side of reaction
            for i in range(2):
                cur_side = sides[i]
                cur_elems = re.findall(r"[A-Z][a-z]?", cur_side)
                # If this is left side, record order of elems
                if i == 0:
                    elems = list(dict.fromkeys((filter(None, cur_elems))))
                    cur_elems = set(elems)
                else:
                    cur_elems = set(filter(None, cur_elems))
                side_elems.append(cur_elems)
                sides[i] = list(filter(None, re.split(r"\s|\+", cur_side)))
            # If left and right side of reaction have matching elems
            if side_elems[0] == side_elems[1]:
                reaction = Reaction(reaction, direction, sides[0], sides[1], elems)
                break;
        error_msg("Invalid input")

    return reaction
                 

# Build coefficient matrix for user chem equation
def build_coeff_matrix(reaction):

    # Initialize vars
    elems = reaction.get_elems()
    num_elems = len(elems)
    num_left_compounds = len(reaction.get_left())
    num_compounds = len(reaction.get_reaction()) - 1
    coeff_matrix = np.zeros((num_elems, num_compounds))
    compound_count = 0;
    cur_elems = []
    elem = ""
    elem_info = []
    sub_coeff = 0;
    row = 0;

    # Iterate through compounds in reaction
    for compound in reaction.get_reaction():
        # If current compound is an actual compound
        if (compound != SEPERATOR):
            cur_elems = re.findall(r"[A-Z][a-z]?[\d]*", compound)
            # Iterate through elems in current compound
            for elem in cur_elems:
                elem_info = list(filter(None, re.split(r"(\d+)", elem)))
                # Find corresponding coeff
                if len(elem_info) == 1:
                    sub_coeff = 1;
                else:
                    sub_coeff = int(elem_info[1])
                # If elem on right side, multiply coeff by -1
                if compound_count >= num_left_compounds:
                    sub_coeff *= -1
                # Update coeff in coeff matrix
                row = elems.index(elem_info[0])
                coeff_matrix[row, compound_count] += sub_coeff
            compound_count += 1
        

    # Iterate first half of equation
    return coeff_matrix


# Solve coeff matrix for coeffs, return fraction and whole number lists
def solve_coeff(coeff_matrix):

    # Initialize vars
    coeff = 0
    cur_frac = 0
    base_coeffs_frac = []
    temp_coeffs = []
    coeffs = []
    
    try:
        # Find null space of coeff matrix, basis serves as basic coeffs
        base_coeffs = linalg.null_space(coeff_matrix)
        print(base_coeffs)
        # Convert float coeffs to fractions
        for coeff in base_coeffs:
            cur_frac = Fraction(float(coeff)).limit_denominator()
            base_coeffs_frac.append(cur_frac)
        # Attempt to remove fractions from coeff matrix
        min_val = base_coeffs.min()
        temp_coeffs = ((1 / min_val) * base_coeffs)
        for coeff in temp_coeffs:
            coeffs.append(coeff[0])
    except:
        error_msg("Equation could not be balanced")
        input(">>> Press <Enter> to exit program ")
        sys.exit()
        
    return (base_coeffs_frac, coeffs)


# Print resultant chem equation with proper coefficients
def print_eq(reaction, coeffs):

    # Initialize vars
    reaction_sections = reaction.get_reaction()
    coeff_indx = 0
    
    # Iterate through compounds in reaction
    for compound in reaction_sections:
        # If current compound is an acutal compound
        if compound != SEPERATOR:
            print(str("({})".format(coeffs[coeff_indx])), end ="")
            coeff_indx += 1
        print(compound, end = " ")
    print()


# Return whether user would like to balance more equations
def more_eq():

    response = input("\nWould you like to balance another " \
                     "equation? (y / n): ").lower()
    result = True if response == "y" else False
    return result


# Print simple invalid input error
def error_msg(msg):

    print(">>> Error: {}".format(msg))


# Equation class storing info on reaction
class Reaction:

    # Equation constructor
    def __init__(self, reaction, direction, left, right, elems):
        self.direction = direction
        self.reaction = reaction
        self.left = left
        self.right = right
        self.elems = elems

    # Return raw reaction data
    def get_reaction(self):
        return self.reaction

    # Return direction of reaction
    def get_direction(self):
        return self.direction

    # Return left side of reaction
    def get_left(self):
        return self.left

    # Return right side of reaction
    def get_right(self):
        return self.right

    # Return elems involved in reaction
    def get_elems(self):
        return self.elems

    # String representation of reaction
    def __str__(self):
        format_str = "Reaction: {}\nDirection: {}\nLeft: {}\nRight: {}\nElems: {}"
        return format_str.format(self.reaction, self.direction,
                                 self.left, self.right, self.elems)


# Run main function
main()
