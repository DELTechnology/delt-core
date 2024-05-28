HOW TO FILL IN EXCEL FILE


Important remarks:
- If your experiment includes the hybridization of two libraries, you have to fill in two separate files (i.e., one file per library).
- Do only edit the sheets "step1", ..., "stepN", and "const".
- Delete the sheets you don't need (e.g., delete "step2" if your experiment includes only one reaction step).
- The sheets "scaffolds" and "smarts" can be ignored.
- Do only use the options available (see below).
- In case of a two step reaction (e.g., SR followed by ABF), write all reaction types separated by a comma (see example usage of step 1).
- If no reaction takes place, leave the corresponding cell empty (see example usage of step 2).



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



#########
# const #
#########

ID
Description:	Number of building block (should match the number of reaction steps)
Options:	1, ..., n

Sequence
Description:	Sequence of the constant regions surrounding the codon of the corresponding building block
Options:	[ATCG]

Reverse
Description:	Whether the final product contains the reverse sequence of the constant region provided
Option:		0, 1

Complement
Description:	Whether the final product contains the complement sequence of the constant region provided
Options:	0, 1

Special
Description:	not yet available
Options:	0


Example usage:

ID	Sequence				Reverse	Complement	Special
-------------------------------------------------------------------------------
1	CTGTGTGCTG{codon}CGAGTCCCATGGCGC	0	1		0
2	CGGATCGACG{codon}GCGTCAGGCAGC		1	0		0



