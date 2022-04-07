# asflux

Python wrapper to run [flux simulator](https://confluence.sammeth.net/display/SIM/Home).

The goals of this project are:
* improve flux simulator userfriendliness
* add AS/LSV simulation (similarly to [AsimulatoR](https://github.com/biomedbigdata/ASimulatoR))
* add support for replicates

The prerequisites are [biopython](https://biopython.org/) and flux-simulator (binary in `$PATH` variable).

```
cd example
tar xvfz data.tar.gz
python3 ../asflux.py 21.fa 21.pc.small.gtf --wd out -l 100 -r 100 
```

### TODO
- [ ] AS events simulation
- [ ] LSV simulation
- [ ] replicates

Nice things to add:
- [ ] allows user to define custom path for flux-simulator
- [ ] more arguments (e.g., polyA, single-end, other flux parameters to tweak simulation)
- [ ] single-end reads support
- [ ] simulate from given .pro
