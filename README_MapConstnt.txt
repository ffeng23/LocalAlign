Notes about mapping constants.

05/27/2018 by reading the code, we try to write notes about the maping_constants

In this program, we use forward and reverse to map constant and UPM with tso. Then we will have different results for mapping, map both, map forward, map reverse, map crossover and map break outside constant.

For the last two results, we don't expect map crossover. This is supposed to be a bad thing. map break outside constant, though, is a good thing. It is a subset of both map both and map forward (here means constant sequences, because 454 in this case are reversed somehow). Most of map both sequences are breakoutside constant, only the case where constant and UPm are connected to each other will be written to map both. It is also possible to have map forward only sequences to fall in break outside constant case. When writting to file, we separate the break outside constant sequences with other cases. There is no overlapping between different files.

***Added comments:
	 It could be the case in terms of map both in which the break could end up with inside the constant region (cross over??). Not necessarilly fall out of constant or exactly end to end connecting to each other with UPM (like mentioned in above). 
