#!/usr/bin/env python
# coding: utf-8

# # Find Unique
#  - Deliverables:
#      - findUnique.py - 50 points ( We will test this code)
#      - Notebook with Inspection materials and Results ( 5 pts )
#      - required input from STDIN
#      - required output to STDOUT
#      - output must be sorted by tRNA header
# 
# ## Sets and tRNA
# 
# Mitochondrial tRNA are encoded in the mitochondrial genome and account for ~10% of the coding space, yet over 50% of mitochondrial genomic disease have their origin in this molecule. These molecules are transcribed, processed, modified, and ultimately folded to produce these mature adapters between the mitochondrial transcription and translation processes. Defects in mt.tRNA maturation reflect additional disease states .. if only we could identify those mutations. We are working on such a device, though we first need to use it to identify abundance of these molecules in a mixed population. Are there any unique subsequences among the 22 mt.tRNA that can be used as "tags"? If so, we can count those tags to assess abundance.
# 
# ## Assignment
# 
# Write a python command line program that reads a file of fasta sequences from STDIN, finds the unique subsequences that occur in each single tRNA such that no members of this set occur among any of the other tRNA sets. Each of those 22 sets should be minimized, such that no member of that unique subsequence set is a substring of any other member of this set. [ Unique and Essential]
# 
# As an example, let's say that both ACG and AAACGA are in a unique set. Since ACG is a substring of AAACGA we would remove AAACGA. [ ACG is Essential ]
# 
# Use Python sets for this assignment[__required__]. Not only will your code be smaller, but it will be more likely to work. The union, intersection and difference operators will be quite useful.
# 
# ## Rough design plan...
# 
# 1) compute the set of all substrings from each tRNA sequence. [powerset] I will refer to a set of substrings as a tRNAset.
# 
# 2) for each tRNAset, compute the union of all other tRNA sets, and then remove that union from the current tRNAset. Notice that this union operation finds all of the substrings from all other tRNA. If any of those are present in your current tRNA, then they are not unique ! [ Unique]
# 
# 3) for each unique tRNAset, it now contains the truly unique ones, along with any extensions of that subsequence. If, for example, it was found that G only occurred in a single tRNA, then adding an A onto that G must also be unique because it has a G in it. We only want the minimal form.. G. [ Essential]
# 
# 4) Remove spaces from the header line before sorting and printing. This will make your output a little prettier.
# 
# 5) make sure to remove any alignment characters from your initial sequences. These characters are periods(.), underscores(\_), or dashes(-) . You will find many new characters in the sequence other than {ACGU} - leave these in place, they are modified bases.
# 
# ## Report
# Print a report that contains items as follows. 
# 
#  - Line 1: the tRNA name
#  - Line 2: the tRNA sequence ( with the alignment characters removed)
#  - lines 3-80 or so, each unique element.
# 
# These unique elements need to be ordered by their starting position in the tRNA sequence, and should be spaced over so they are directly under that subsequence in the original tRNA. This looks like an alignment, but you can find where it belongs by using the string.find() method. Include dots in all positions to the left, serving as alignment guides for your reader. [ see sample output below ]
# 
# Do this for all of the 22 tRNA sequences.
# 
# Print the tRNA out as above, sorted by the header line.  
# 
# ## Hints:
# 
# __use sets and the set operators - this is required !__
# 
# your final code will be under 100 lines.
# 
# Do most of your coding using methods defined in a class and then instantiated as objects. 
# 
# The sequences include characters that are just alignment characters. They are not part of the sequence and must be removed. [-\_\.] are alignment characters. 
# 
# When removing items from a tRNAset, don't do this while iterating through that set. Also, when building a unique set you will need the original contents of all other tRNAsets. So.. build a new set to keep the unique contents, or keep track of the elements you intend to delete later. Notice that you build the union of all other tRNA, and this happens 22 times - these unions are all distinct from each other. Example, consider 4 sets, A,B,C,D.  we would compute the set B, C, D to use against set A, and we would compute the set A, C, D to use against B. A and B need to not change while we are doing this computation.
# 
# There are cases where a unique element is present multiple times in a single tRNA. This sequence is unique to this tRNA and should be reported in its multiple locations ( same sequence on multiple lines of the output).
# 
# This is a command line program, though it does not have any optional parameters. You really don't need the commandline class. Your input comes from STDIN and output goes to STDOUT. Use the FastaReader for input and print statement for output. If you are going to use jupyter to develop your code, you can add a filename to fastareader for your testing.
# 
# ## Submission
# 
# Submit code using Canvas. As always, work with your group and write your own code.

# In[2]:


#!/usr/bin/env python3
# Name: Yashesha Kothari (ykothari)
# Group Members: Srikar Bevara, Cameron Ahayee, Donya Mirzadeh, Tia Abraham

import sys

class FastAreader:
    def __init__(self, fname=''):
        '''constructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta(self):
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

            yield header, sequence

class FindUnique:
    """This class reads a FASTA file of tRNA sequences and finds the unique subsequences for each tRNA using power set method.
    Removes all duplicate subsequences found in all other tRNAs except the one current.
    """
    def __init__(self):
        """initializes the instance variables for FindUnique"""
        
        self.powerSetList = [] #list with all power sets for each tRNA
        self.uniqueList = [] #list with unique subsequences for each tRNA
        self.tRNADict = {} #dictionary containing the header and seqeunces for each tRNA
        self.fasta_file = FastAreader('bos-tRNA.fa') #object to read in the FASTA file
        count = 0
        
        for header, sequence in self.fasta_file.readFasta():
            #iterate through header and sequence of FASTA file 
            
            newSequence = self.deleteChar(sequence) #uses deleteChar method and stores
            self.tRNADict[count] = [header, newSequence] #creates an entry in dictionary for current count with header and newSquence
            myPowerSet = self.powerSet(newSequence) #generates power set 
            self.powerSetList.append(myPowerSet) #appends power set to List
            count += 1
            
    def deleteChar(self, sequence):
        """For each sequence, removes dashes and underscores"""
        
        no_dash_under = sequence.replace('_','').replace('-','')
        return no_dash_under
    
    def powerSet(self, sequence):
        """Returns the power set of a given sequence"""
        powerSet = set() #creates an empty set to store the power set
        
        #Loop through the input sequence
        for i in range(len(sequence)):
            length = len(sequence) #set the length of sequence
            
            # While the length is greater than i, add the subsequence from index i to length to the power set.
            while length > i:
                powerSet.add(sequence[i : length])
                length -= 1
        return powerSet
    
    def findUniques(self):
        """
        Finds the unique tRNA subsequences by removing duplicate sequences from all but the appropriate set.
    Modifies the 'uniqueList' attribute of the object.
        """
        
        #Loop through each power set in the list attribute
        for copySet in self.powerSetList:
            union = set() #create an empty set to store union of all the other tRNA sets
            
            #Make a copy of the list attribute and remove current power set 
            copyList = self.powerSetList.copy()
            copyList.remove(copySet)  # Remove power set from list.
            
            #loop through each power in copy list and add to union
            for powerSet in copyList:
                union = union.union(powerSet)  # The union of all the other tRNA sets.
            copySet.difference_update(union)  # Update the copySet, removing elements found in union.

            #create new set as the copy
            newSet = copySet.copy()
            
            #loop thorugh each string in current power set
            for string1 in copySet:
                
                #create copy of current power set and remove current string
                uniqueSet = copySet.copy()
                uniqueSet.remove(string1)
                
                #loop through each string in unique set and check if current string is substring
                for string2 in uniqueSet:
                    if string1 in string2 and len(string1) < len(string2):
                        newSet.discard(string2)  # remove larger

            self.uniqueList.append(newSet) #add
        
    def printSequences(self):
         """Returns the power set of a given sequence"""
        
        #loop thorugh each tRNA sequence in the 'tRNADict' attribute
        for i in range(0,len(self.tRNADict)):
            
            #Get the header and sequence of current tRNA sequence
            mySequence = self.tRNADict[i]
            header = mySequence[0]
            sequence = mySequence[1]
            
            #print each
            print(header)
            print(sequence)
            
            #Get length and loop though each unique substring in current tRNA
            length = len(sequence)
            for place in range(0,length):
                for substring in self.uniqueList[i]:
                    substringLength = len(substring)
                    
                    #if current substring matches the substring in the sequence starting from current position
                    if substring == sequence[place : place + substringLength]:
                        result = ('.') * place + substring #create output with dots for positions before
                        print (result)

def main():
    """
    Runs FindUnique class to identify unique subsequences in tRNA sequences.
    Print the sequences.
    """
    out_tRNA = FindUnique()
    out_tRNA.findUniques()
    out_tRNA.printSequences()

if __name__ == "__main__":
    main()          
        


# # Sample Output

# In[ ]:


tRNA|Lys|∃UU|Bostaurus|mitochondrial
CACUAAGA"LCUAUAUAGCACPAACCU∃UU6AGUUAGAGAUUGAGAGCCAU"UACUCUCCUUGGUGACCA
CACU
.ACUA
...UAA
....AAG
.......A"
.........L
..........CUAU
............AUAU
..............AUAG
...............UAGC
.................GCA
.....................P
......................AAC
.......................ACCU
...........................∃
..............................6
...............................AGU
................................GUU
.................................UUA
..................................UAGA
...................................AGAGA
......................................GAU
.......................................AUU
........................................UUGA
.........................................UGAG
..........................................GAGAG
............................................GAGC
..............................................GCC
................................................CAU
..................................................U"
..................................................."U
....................................................UAC
.....................................................ACUC
.......................................................UCU
.........................................................UCC
...........................................................CUU
..............................................................GG
...............................................................GUG
.................................................................GAC
..................................................................ACCA


# My inspection team is Srikar, Donya, Cameron, and Tia. They found:
# 
# - I had a lot of repetitive unecessary lines of code that could have been narrowed down, so I combined some and deleted some altogether
# - did not put the file name as the FastAreader attribute so it was not running
# - the deleteChar function could have been implemented into the __init__ but i decided to leave it there because easier for me to follow
# - there are many other design options that could make this more straightforward
# - I did not have docstrings for the functions and could add more specific comments
# 
# I fixed most of the suggestions and added more comments and docstrings

# In[ ]:




