f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/code1_NF2.txt'
with open(f, 'r') as f:
    c1 = f.read().split('\n')
c1 = set(c1)

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/code2_NF2.txt'
with open(f, 'r') as f:
    c2 = f.read().split('\n')
c2 = set(c2)

c1.intersection(c2)

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/NF2_primers1b.txt'
with open(f, 'r') as f:
    c3 = f.read().split('\n')
c3 = set(c3)

f = '/Users/adrianomartinelli/projects/delt/DECL-data-analysis/simulation/PrimerCodes1.txt'
with open(f, 'r') as f:
    c4 = f.read().split('\n')
c4 = set(c4)

# not true, i.e. we have duplicates
len(c1.union(c2).union(c3).union(c4)) == len(c1) + len(c2) + len(c3) + len(c4)

# here we have an intesection
c1.intersection(c4)
