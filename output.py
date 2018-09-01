"""
Stuff to do with emitting output. Pretty-printers and so on. Things not to do
with actually solving the problem.
"""

from termcolor import colored

def printErrors(errors, read, ref, colouredDiff):
    """
    Output the errors for this read. Either as a coloured diff, or in a text summary.

    @param errors A list of error objects.
    @param read The read, post-alignment, a Seq object.
    @param ref The reference, post-alignment, a Seq object..
    @param colouredDiff If true, will print a fancy coloured diff instead of a concise error summary.
    """

    print(errors)
    
    if colouredDiff:
        print()
        # Print a coloured diff :D
        for i in range(0, 100, 5):
            print(str(i) + "   ", end='')
            if len(str(i)) == 1:
                print(" ", end='')

        print()
        print("|...." * 20)

        for i in range(0, len(ref)):
            if i > 0 and i % 100 == 0:
                print()  # Add a newline after each 100 symbols.

            refChar = str(ref)[i]
            readChar = str(read)[i]

            if refChar != "-":
                if refChar == readChar:
                    print(colored(refChar, 'green'), end='')  # no change
                elif readChar != "-":
                    print(colored(readChar, 'red'), end='')  # Substitution
                else:
                    print(refChar, end='')  # Deletion
            else:
                print(colored(readChar, 'blue'), end='')  # Insertion
                
        print()

def print_coloured_diff(errors, read, ref, PRINT_COLOURED_DIFF):
    print("\n\n#############################################################################")
    printErrors(errors, read, ref, PRINT_COLOURED_DIFF)

def printSummary(dictionary, keys):
    print(','.join(keys))
    for mutation, counts in sorted(dictionary.items()):
        s = [str(counts.get(library)) for library in keys[1:]]
        print(mutation, ','.join(s), sep=",")