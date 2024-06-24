
Defining an input model
#######################

Consider a toy gene regulatory network with a series of mutual
inhibition motifs as illustrated below:

.. figure:: figs/boolode-example.png
   :align: center

   If there is an equal probabiliy of taking either the left hand branch
   or the right hand branch starting at (g1=high,g2=high), this network will
   result in 8 distinct steady states, with one of the eight extreme genes being
   active in each steady state.

BoolODE can be used to generate gene expression data of a sample of cells undergoing
a developmental process controlled by a gene regulatory network with a wiring diagram
as shown above.

The above network diagram can be represented as a Boolean model as follows:
   
.. code:: text
          
          Gene	Rule
          g1	( g1 ) and not ( g2 )
          g2	( g2 ) and not ( g1 )
          g11	( g1 or g11 ) and not ( g12 )
          g12	( g1 or g12 ) and not ( g11 )
          g21	( g2 or g21 ) and not ( g22 )
          g22	( g2 or g22 ) and not ( g21 )
          g111	( g11 or g111 ) and not ( g112 )
          g112	( g11 or g112 ) and not ( g111 )
          g121	( g12 or g121 ) and not ( g122 )
          g122	( g12 or g122 ) and not ( g121 )
          g211	( g21 or g211 ) and not ( g212 )
          g212	( g21 or g212 ) and not ( g211 )
          g221	( g22 or g221 ) and not ( g222 )
          g222	( g22 or g222 ) and not ( g221 )


In order to test BoolODE, copy the rules  above into a file, say ``multistate.txt``,
in the ``/data`` folder in the BoolODE folder.

.. note:: BoolODE expects a tab separated rule file! Make sure the
           input file has tabs between the ``Gene`` and ``Rule`` column.
          
Notice that the ``g1<->g2`` interaction is central to the behavior of
the model.  Let us choose an initial condition such that these genes
have a medium/high value. Biologically, this can be interpreted as the
initial undifferentiated cells have genes ``g1`` and ``g2`` 'ON', and
all other genes 'OFF'.  We can specify this initial state by creating a file with the
following lines:

.. code:: text

          Genes	Values
          ['g1','g2']	[1,1]

The first column takes a list of genes formatted like a python list of
strings. The second column takes a list of values to assign to the
list of genes, again formatted as a python list of values. **Any gene not
listed in this file will be assigned an initial value of 0 by default.**

Let us save the initial condition specification in
``/data/multistate_ics.txt``.



.. note:: The initial condition file should also be tab
           separated. Also, pay attention to the capitalization
           of the column names.


.. tip:: We recommend starting off with initial condition values
         between 0 and 2, since larger values might cause numerical
         instability.


Running BoolODE
###############

A few more decisions before we get our simulated dataset of gene expression
values:

1. How long do you want to simulate this model? Let's choose 8 time units
2. How many cells do you want to simulate? 1000 seems like a good number.
3. Where would you like to store the simulation output? Say, ``./test/``. BoolODE will create any folders that don't already exist.   
4. Do you want to speed up the simulations by running simulations in parallel? Yes!
5. Finally, we can process the output from BoolODE and generate a sample of single cells
   by specifying.

Users can specify simulation settings using a YAML file. Such a file should have
three sections:
1. ``global_settings`` defines the paths that BoolODE needs to read inputs and write outputs, as well as model type to generate
2. ``jobs`` specifies the list of models to simulate. Right now we have one model, for which we have defined the parameters for simulation.
3. ``post_processing``  specifies how BoolODE will process the output simulations in order to serve as input files to BEELINE

A sample config file with the settings above would look something like the following.

.. code:: yaml
          
          global_settings:                                  
            model_dir: "data"                               
            output_dir: "test"                              
            do_simulations: True                       
            do_post_processing: False
            modeltype: 'hill'                               
          jobs:                                             
            - name: "multistate-1"                          
              model_definition: "multistate.txt"            
              model_initial_conditions: "multistate_ics.txt"
              num_cells: 500                                
              simulation_time: 8                            
              do_parallel: True                             
          post_processing:                                  
            genSamples:                                     
              - sample_size: 400                            
                nDataSets: 1                                
            DimRed:                                         
              - perplexity: 100                             
          
 For a full list of available options, see `example-config.yaml`

Working with BoolODE output
###########################

The simulations will take about 5 minutes to complete. At the end of a successful run, the
output directory should like this:

.. code:: text

          test/                      # User specified destination
          |-- parameters.txt         # Parameter names and values generated for input model
          |-- ExpressionData.csv     # Gene expression dataset
          |-- PseudoTime.csv         # Simulation time of each sample time point/cell
          |-- refNetwork.csv         # Boolean network represented as an edge list, the ground truth network
          `-- simulations/
              |-- E1.csv
              |-- E2.csv
              ...

Where ``E1.csv, E2.csv, ...`` are individual simulations. Each column in these
files, the cell IDs, has the form ``E<simulation number>_<timepoint>``.
          
The above run has generated a series of simulation files. In the next step, we
can use the same config file to carry out post processing on the simulation output,
as described in :ref:`geninputs`.

..
   Notice the eight steady state clusters! This dataset can now be
   processed further using the tools described in  :ref:`geninputs` to produce
   input datasets for the BEELINE pipeline.
            
