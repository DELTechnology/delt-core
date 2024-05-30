HOW TO FILL IN EXCEL FILE


Important remarks:
- If your experiment includes the hybridization of two libraries, you have to fill in two separate Excel files (i.e., one file per library).
- Do only edit the sheets "step1", ..., "stepN", "scaffolds", and "const" (i.e., the sheet "smarts" can be ignored).
- Delete the sheets you don't need (e.g., delete "step2" if your experiment includes only one reaction step).
- Do only use the options available (see below).
- In case of a two step reaction (e.g., SR followed by ABF), write all reaction types separated by a comma (see example usage of step 1).
- If no reaction takes place or no scaffold is used, leave the corresponding cell empty (see example usage of step 2).
- Only reactions between two building blocks or between a building block and a scaffold are considered (i.e., reactions between a building block and a DNA strand should be neglected).
- All sequences provided should be written in the 5'-to-3' direction.
- In case of a hybridization, the constant regions provided must not overlap (see example usage below).



##########
# STEP 1 #
##########

ID
Description:	ID of the first building block
Options:	1, ..., n

SMILES
Description:	SMILES of building block n
Options:	Valid SMILES structure

ScaffoldID
Description:	Iodine position (ortho, meta, para)
Options		o, m, p

Codon
Description:	Codon of building block n
Option:		[ATCG]

ReactionType
Description:	Reaction type(s) (amide bond formation, Staudinger reduction, cycloaddition, Suzuki, Sonogashira, dehalogenation)
Options:	ABF, SR, CuAAC, Suz, Son, DH


Example usage:

ID	SMILES		ScaffoldID	Codon	ReactionType
------------------------------------------------------------
1	NCC#C		p		GTGGTG	CuAAC
2	OCC(O)=O	m		ATATTG	SR,ABF
3	N		o		TAAGTA	SR



##################
# STEP 2, ..., N #
##################

ID
Description:	ID of the N-th building block
Options:	1, ..., n

SMILES
Description:	SMILES of building block n
Options:	Valid SMILES structure

Codon
Description:	Codon of building block n
Option:		[ATCG]

ReactionType
Description:	Reaction type(s) (amide bond formation, Staudinger reduction, cycloaddition, Suzuki, Sonogashira, dehalogenation)
Options:	ABF, SR, CuAAC, Suz, Son, DH


Example usage:

ID	SMILES		Codon	ReactionType
--------------------------------------------
1	OB(O)c1occc1	TGCCTTC	Suz
2	C#CC1CCCCC1	TGGTAGA	Son
3	[I]		TGGCAAT



#############
# SCAFFOLDS #
#############

ScaffoldID
Description:	ID of the corresponding scaffold
Options:	1, ..., n

SMILES
Description:	SMILES of scaffold n
Options:	Valid SMILES structure


Example usage:

ScaffoldID	SMILES
------------------------------------------------
1		Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1
2		[N-]=[N+]=NC(C(O)=O)Cc1cc(I)ccc1
3		[N-]=[N+]=NC(C(O)=O)Cc1c(I)cccc1



#########
# CONST #
#########

Sequence
Description:	Sequence of the constant regions surrounding the codons of the building blocks
Options:	[ATCG]

Reverse
Description:	Whether the final DNA sequence contains the reverse sequence of the constant region provided
Option:		0, 1

Complement
Description:	Whether the final DNA sequence contains the complement sequence of the constant region provided
Options:	0, 1


Example usage (one library):

Final DNA sequence:	5' GAC{codon}GCT{codon}GCT 3'

Sequence			Reverse	Complement
--------------------------------------------------
GAC{codon}GCT{codon}GCT		0	0


Example usage (hybridization of two libraries):

5' GTA{codon}GGA              3'
   |||       |||
3' CAT{xxxxx}CCTGAA{codon}TGC 5'

Final DNA sequence:	5' GTA{codon}GGACTT{codon}ACG 3'

Excel file #1:
Sequence				Reverse	Complement
----------------------------------------------------------
GTA{codon}GGA				0	0

Excel file #2:
Sequence				Reverse	Complement
----------------------------------------------------------
CGT{codon}AAG				1	1

