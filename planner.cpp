/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <cmath>
#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 
#include <unordered_set>
#include <limits>

const int INF = std::numeric_limits<int>::max(); // Infinity

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTSTAR     1
#define PRM         2

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                GIVEN FUNCTIONS                                                    //
//                                                                                                                   //
//*******************************************************************************************************************//

/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                          DEFAULT PLANNER FUNCTION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
		
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf) {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;
    
    return;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              RRT IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRT(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */

    planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}


//*******************************************************************************************************************//
//                                                                                                                   //
//                                           RRT STAR IMPLEMENTATION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void plannerRRTStar(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
    planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              PRM IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

double distanceCalc(double *anglesV_rad1,double *anglesV_rad2,int numOfDOFs){
	double distance=0;
    for (int j = 0; j < numOfDOFs; j++){
        if(distance < fabs(anglesV_rad1[j] - anglesV_rad2[j]))
            distance = fabs(anglesV_rad1[j] - anglesV_rad2[j]);
    }
	return distance;
}

static void plannerPRM(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
	const int N = 100; // Number of random configurations to generate
	const int K=10; //nearest neighbours
	//no plan by default
	*plan = NULL;
	*planlength = 0;
    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, PI); // Range for degrees of freedom: [0, π]

    double V[N+2][numofDOFs]; //set of vertices
    int count = 0;
    while (count < N) {
        // Generate random configuration
        double random_config[numofDOFs];
        for (int i = 0; i < numofDOFs; ++i) { 
            random_config[i] = round(dis(gen) * 1000) / 1000; //this gives output of 3 decimal places
        }
        // Check if the configuration is valid
        if (IsValidArmConfiguration(random_config, numofDOFs, map, x_size, y_size)) {
            for (int i = 0; i < numofDOFs; ++i) {
                V[count][i] = random_config[i];
            }
		}
            count++;
    }
	double distance_matrix[N][N];
	// declaring the distacne matrix to have all elements -1
	for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            distance_matrix[i][j] = -1;
        }
    }
	//calculating all distances between i and j (halfing it due to it being a symmetrical matrix)
	
	for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
			if(distance_matrix[j][i]!=-1){
				distance_matrix[i][j] = distance_matrix[j][i];
			}
			else{
            distance_matrix[i][j] = distanceCalc(V[i], V[j],numofDOFs);
			}
        }
    }

	int nearest_neighbors[N+2][N+2];
    for (int i = 0; i < N; ++i) {
        // Initialize nearest_neighbors matrix with INF
        for (int j = 0; j < N; ++j) {
            nearest_neighbors[i][j] = INF;
        }
	}
	//sorting the elements and then choosing the closest points 
    for (int i = 0; i < N; ++i) {
        // Create an array to store the indices of the points sorted by distance
        int sorted_indices[N];
        for (int j = 0; j < N; ++j) {
            sorted_indices[j] = j;
        }

        // Sort the indices based on the distance from the current point
        std::sort(sorted_indices, sorted_indices + N, [&](int a, int b) {
            return distance_matrix[i][a] < distance_matrix[i][b];
        });

        // Store the indices of the K nearest neighbors
        for (int k = 0; k < K; ++k) {
            nearest_neighbors[i][sorted_indices[k]]=distance_matrix[i][sorted_indices[k]];
        }
    }
	int maxsizeedge=N*K;
	double E[maxsizeedge][2];
	//this nearest_neighbours matrix is already an adjaceny matrix for the graph 
	//(but all lines here are considered eligible so we can now check if the lines paths on these edges are eligible
	int countedge=0;
	for (int i = 0; i < N; ++i) {
        for(int k=0;k <N;++k){
			if(nearest_neighbors[i][k]==INF){
				continue;
			}
			//coordinates of base for final and initial positoions of the motion
			double xi= ((double)x_size)/2.0;
			double yi = 0;
			double xf=((double)x_size)/2.0;
			double yf=0;
			for(int j = 0; j < numofDOFs; j++){
			//compute the corresponding line segment for end point of each link 
				xi = xi + LINKLENGTH_CELLS*cos(2*PI-V[i][j]);
				yi = yi - LINKLENGTH_CELLS*sin(2*PI-V[i][j]);
				xf = xi + LINKLENGTH_CELLS*cos(2*PI-V[k][j]);
				yf = yi - LINKLENGTH_CELLS*sin(2*PI-V[k][j]);
				//check the validity of the corresponding line segment
				if(IsValidLineSegment(xi,yi,xf,yf,map,x_size,y_size)){
					E[countedge][0]=i;
					E[countedge][1]=k;	
				}			
			}    
		}
	}
	for(int i=0;i<numofDOFs;i++){
		V[N][i]=armstart_anglesV_rad[i]; //adding the start point
		V[N+1][i]=armgoal_anglesV_rad[i]; //adding the end point 
	}
	//we need to now fill the neareast neighbour matrix for the start and end point
	for(int i=0;i<(N+2);i++){
		nearest_neighbors[i][i]=INF;
	}
	double maxdist=0;
	//calculating distance
	for(int i=N;i<=(N+1);i++){
		for(int j=0;j<(N+2);j++){
			nearest_neighbors[i][j]=distanceCalc(V[i], V[j],numofDOFs);			
		}
	}
	for(int i=N;i<=(N+1);i++){
		for(int j=0;j<(N);j++){
			double maxj=-1;
			for(int k=0;k<N;k++){
				if (nearest_neighbors[j][k]>maxj && nearest_neighbors[j][k]!=INF){
					maxj=nearest_neighbors[j][k];
				}
			}
			if(maxj==-1 || maxj==INF){
				nearest_neighbors[i][j]=INF;
			}
		}
	}	

	// int sorted_indices[N+2];
    // for (int j = 0; j < (N+2); ++j) {
    //     sorted_indices[j] = j;
    // }
	// for(int i=N;i<=(N+1);i++){
	// 	std::sort(sorted_indices, sorted_indices + N, [&](int a, int b) {
	// 		return nearest_neighbors[i][a] < nearest_neighbors[i][b];
	// 	});
    //     for (int k = K; k < N+2; ++k) {
    //         nearest_neighbors[i][sorted_indices[k]]=INF;
    //     }
	// }
// Function to implement Dijkstra's algorithm
	int distanceArr[(N+2)];
	bool visited[N+2];
	int parent[N+2];

	// Initialize distances, visited array, and parent array
	for (int i = 0; i < (N+2); ++i) {
		distanceArr[i] = INF;
		visited[i] = false;
		parent[i] = -1;
	}
	distanceArr[N+1] = 0;

	// Main loop
	for (int count = 0; count < (N+2) - 1; ++count) {
		int minDist = INF, minIndex = -1;
		for (int v = 0; v < (N+2); ++v) {
			if (!visited[v] && distanceArr[v] < minDist) {
				minDist = distanceArr[v];
				minIndex = v;
			}
		}
		int u = minIndex;
		if (u == -1) break;
		visited[u] = true;
		// Update distances to neighbors of u
		for (int v = 0; v < (N+2); ++v) {
			if (!visited[v] && nearest_neighbors[u][v] && distanceArr[u] != INF && distanceArr[u] + nearest_neighbors[u][v] < distanceArr[v]) {
				distanceArr[v] = distanceArr[u] + nearest_neighbors[u][v];
				parent[v] = u;
			}
		}
	}
	int revpath[N+2];
	int current = N+1;
	int index = 0;
	int pathcount=0;
	while (current != -1) {
		revpath[index++] = current;
		current = parent[current];
		pathcount++;
	}
	if (current!=N){
		printf("ERROR: no path");
	}
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i <pathcount; i--){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = V[revpath[pathcount-i-1]][j];
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf) {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = pathcount;
    
    /* TODO: Replace with your implementation */

    //planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
	return;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                MAIN FUNCTION                                                      //
//                                                                                                                   //
//*******************************************************************************************************************//

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
		printf("Invalid start or goal configuration");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;

    // Call the corresponding planner function
    if (whichPlanner == PRM)
    {
        plannerPRM(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRT)
    {
        plannerRRT(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRTSTAR)
    {
        plannerRRTStar(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else
    {
        planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
	printf("We did it -part1");
	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my 
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////
	

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
		printf("Start or goal postion not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
