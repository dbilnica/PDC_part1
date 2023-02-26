Class city {
Int index;
double lb; // two mins, select lower one
Path paths[] = {Path2, Path3…}
}

Class Path {
City first;
City second;
double cost;
}

tour = {0,2,3,1} - indexes of visited cities, with tour we also save current tours cost

For path in paths:
	Calculate lower bound for each path of the city
Then move to next city

* Use double for everything
* Every city can be visited only once
* Tour, ukladam si indexy kde jsem byl (budu vedet, ze tama uz nemam jit)
* Prunning should be in the queue
* 