

  -h, --help            show this help message and exit
  --max-time=MAX_TIME   Total time of simulation. (Default = 20)
  --num-cells=NUM_CELLS
                        Number of cells sample. (Default = 100)
  --sample-cells        Sample a single cell from each trajectory? By default
                        will store full trajectory of each simulation (Default
                        = False)
  --add-dummy           Add dummy genes
  -n, --normalize-trajectory
                        Min-Max normalize genes across all experiments
  -i, --identical-pars  Set single value to similar parameters NOTE: Consider
                        setting this if you are using --sample-pars.
  -s, --sample-pars     Sample rate parameters around the hardcoded means ,
                        using 10% stand. dev.
  --std=STD             If sampling parameters, specify standard deviation.
                        (Default 0.1)
  --write-protein       Write both protein and RNA values to file. Useful for
                        debugging.
  -b, --burn-in         Treats the first 25% of the time course as burnin ,
                        samples from the rest.
  --outPrefix=OUTPREFIX
                        Prefix for output files.
  --path=PATH           Path to boolean model file
  --inputs=INPUTS       Path to input parameter files. This is different from
                        specifying a parameter set!
  --pset=PSET           Path to pre-generated parameter set.
  --ics=ICS             Path to list of initial conditions
  --strengths=STRENGTHS
                        Path to list of interaction strengths
  --species-type=SPECIES_TYPE
                        Path to list of molecular species type file.Useful to
                        specify proteins/genes
  -c NCLUSTERS, --nClusters=NCLUSTERS
                        Number of expected clusters in the dataset. (Default =1)
  --max-parents=MAX_PARENTS
                        Number of parents to add to dummy genes. (Default = 1)
  --do-parallel         Run simulations in parallel. Recommended for > 50
                        simulations
