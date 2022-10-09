# asflux

Python wrapper to enhance [flux simulator](https://confluence.sammeth.net/display/SIM/Home).

The goals of this project are:
* improve flux simulator userfriendliness
* add AS/LSV simulation (similarly to [AsimulatoR](https://github.com/biomedbigdata/ASimulatoR))
* add support for replicates

### Installation
```
git clone https://github.com/ldenti/asflux.git
cd asflux
mamba create -c bioconda -n asflux biopython gffutils suppa openjdk=8.0.332 flux-simulator
```

### Example
```
cd example
tar xvfz data.tar.gz
python3 ../asflux.py 21.fa 21.pc.small.gtf -n 10
```

### TODO
- [ ] 2 conditions + replicates
- [ ] allows user to define custom path for suppa.py and flux-simulator
- [ ] more arguments (e.g., polyA, single-end, other flux parameters to tweak simulation)
- [ ] allow chromosomes directory
- [ ] simulate from given .pro
- [ ] improve choice of genes that undergo alternative splicing
- [ ] multiple events per gene
- [ ] LSV simulation
