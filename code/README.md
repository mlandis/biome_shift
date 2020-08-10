To run the RevBayes analysis, call `rb run_inf.Rev` from within this directory. That script will in turn call the `model_biome.Rev` script, which builds the time-stratified regional biome shift model. There are many settings you can tweak within `run_inf.Rev`. I've done my best to comment the scripts so that others may modify them freely. Feel welcome to contact me with any quesitons.

The `supplement_1.Rev` script computes the rate matrix Q(1) described in Supplement 1.