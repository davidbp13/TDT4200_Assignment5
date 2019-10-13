
#include <getopt.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <complex.h>
#include <pthread.h>
#include "libs/bitmap.h"
#include "libs/utilities.h"
#include "mandelCompute.h"
#include "mandelColours.h"

/*
   Producer/Consumer Scheme
   job struct has a function pointer and parameter with which this function is called
   job queue is a linked list, which is empty if the list head is NULL.
   Each element points to the next or NULL if last element
   functions putJob and popJob are inserting and taking a job to/from the jobQueue.
   Those functions are NOT thread safe and must be called only one at a time.
   */

typedef struct job {
	void (*callback)(dwellType *, unsigned int const, unsigned int const, unsigned int const);
	dwellType *dwellBuffer;
	unsigned int atY;
	unsigned int atX;
	unsigned int blockSize;
} job;

typedef struct jobQueue {
	job item;
	struct jobQueue *next;
} jobQueue;


// On Error jobQueue will be freed and the application should exit
void putJob (jobQueue **head, job newJob) {
	jobQueue *new = malloc(sizeof(jobQueue));
	if (new != NULL) {
		new->item = newJob;
		new->next = *head;
		*head = new;
	} else {
		fprintf(stderr,"Could not allocate memory for new job!");
		while ((*head) != NULL) {
			jobQueue *tmp = (*head)->next;
			free(*head);
			(*head) = tmp;
		}
	}
}

// Only allowed to be called if something is in the job queue!
job popJob (jobQueue **head) {
	jobQueue *current = *head;
	job poppedJob = current->item;
	(*head) = current->next;
	free(current);
	return poppedJob;
}

jobQueue *jobQueueHead = NULL;
dwellType *dwellBuffer = NULL;

void createJob(void (*callback)(dwellType *, unsigned int const, unsigned int const, unsigned int const),
			   dwellType *buffer,
			   unsigned int const atY,
			   unsigned int const atX,
			   unsigned int const blockSize)
{
	job newJob = { .callback = callback, .dwellBuffer = buffer, .atY = atY, .atX = atX, .blockSize = blockSize };
	putJob(&jobQueueHead, newJob);
}


/*
   Now the 2 two functions are following which do the computation
   marianiSilver is a subdivsion function computing the Mandelbrot set
   escapeTime is the traditional algorithm which you already know
   */

bool markBorders;
unsigned int blockDim;
unsigned int subdivisions;
bool useMarianiSilver = false;
unsigned int useThreads = 4;

void marianiSilver( dwellType *buffer,
					unsigned int const atY,
					unsigned int const atX,
					unsigned int const blockSize)
{
	dwellType dwell = commonBorder(buffer, atY, atX, blockSize);
	if ( dwell != dwellUncomputed ) {
		fillBlock(buffer, dwell, atY, atX, blockSize);
		if (markBorders)
			markBorder(buffer, dwellBorderFill, atY, atX, blockSize);
	} else if (blockSize <= blockDim) {
		computeBlock(buffer, atY, atX, blockSize);
		if (markBorders)
			markBorder(buffer, dwellBorderCompute, atY, atX, blockSize);
	} else {
		// Subdivision
		unsigned int newBlockSize = blockSize / subdivisions;
		for (unsigned int ydiv = 0; ydiv < subdivisions; ydiv++) {
			for (unsigned int xdiv = 0; xdiv < subdivisions; xdiv++) {
				marianiSilver(buffer, atY + (ydiv * newBlockSize), atX + (xdiv * newBlockSize), newBlockSize);
			}
		}
	}
}

void escapeTime( dwellType *buffer,
				 unsigned int const atY,
				 unsigned int const atX,
				 unsigned int const blockSize)
{
	computeBlock(buffer, atY, atX, blockSize);
	if (markBorders)
		markBorder(buffer, dwellBorderCompute, atY, atX, blockSize);
}

pthread_mutex_t mutex;
pthread_cond_t cond;

void *worker(void *id) {
	long tid;
    tid = (long)id;
    //printf("Hello World! It's me, thread #%ld!\n", tid);
    job threadJob;
    
	// Consume 
	// Take a job from the jobQueue
	while(jobQueueHead != NULL) {
        pthread_mutex_lock(&mutex); //Start by locking the queue mutex
       // printf("It's me, thread #%ld and I got the mutex\n", tid);

        // If the jobQueue is empty, wait for jobs to show up
        //while(jobQueueHead == NULL) {
            // Wait until a new job is in the queue
            //pthread_cond_wait(&cond, &mutex);
        //}

        threadJob = popJob(&jobQueueHead); //As we returned from the call, there must be new data in the queue - get it
        pthread_mutex_unlock(&mutex); //Now unlock the mutex
        
        // Process the job that was taken
        if (useMarianiSilver) {
			/*// Scale the blockSize from res up to a subdividable value
			// Number of possible subdivisions:
			unsigned int const numDiv = ceil(logf((double) resolution/blockDim)/logf((double) subdivisions));
			// Calculate a dividable resolution for the blockSize:
			unsigned int const correctedBlockSize = pow(subdivisions,numDiv) * blockDim;
			// Mariani-Silver subdivision algorithm
			marianiSilver(dwellBuffer, 0, 0, correctedBlockSize);*/
		} 
		else {
			// Traditional Mandelbrot-Set computation or the 'Escape Time' algorithm
			// computeBlock respects the resolution of the image, so we scale the blocks up to
			// a divideable dimension
			//printf("It's me, thread #%ld and I'm pocessing a job. Y = %d, X = %d, size = %d.\n", tid, threadJob.atY, threadJob.atX, threadJob.blockSize);
			threadJob.callback(threadJob.dwellBuffer, threadJob.atY, threadJob.atX, threadJob.blockSize);
			//printf("It's me, thread #%ld and I finished my job\n", tid);
			printf("It's me, thread #%ld and I'm finished a job. Y = %d, X = %d, size = %d.\n", tid, threadJob.atY, threadJob.atX, threadJob.blockSize);
		}
    }    
    printf("It's me, thread #%ld and I exited the while loop\n", tid);
    
    pthread_exit(NULL);
}

void initializeWorkers(unsigned int threadsNumber) {
	pthread_t threads[threadsNumber]; // Array holding each of the threads
	pthread_attr_t attr;			  // Threads attribute
	int rc;							  // Return value for threads creation
	long t;							  // Thread iterator
	
	// Make threads joinable
	pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
	// Initialize the mutex and the condition variable
	pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);
    
    // In the beggining there are no jobs in the queue. An  initial job or jobs are created depending on the algorithm
    if (useMarianiSilver){
		createJob(marianiSilver, dwellBuffer, 0, 0, resolution * resolution); // For Maniani-Silver one job is created and this job will create more as it is necessary
	}
	else{
		// For Escape-Time threadsNumber * threadsNumber jobs are created at once
		unsigned int block = ceil((double) resolution / threadsNumber);
		for (unsigned int t = 0; t < threadsNumber; t++) {
			for (unsigned int x = 0; x < threadsNumber; x++) {
				createJob(escapeTime, dwellBuffer, t * block, x * block, block);
			}
		}	
	}
	
	// Create threads
	for(t=0; t<threadsNumber; t++){
       printf("In main: creating thread %ld\n", t);
       rc = pthread_create(&threads[t], NULL, worker, (void *)t);
       if (rc){
          printf("ERROR; return code from pthread_create() is %d\n", rc);
          exit(-1);
       }
    }
    
    // Wait for threads to complete
    pthread_attr_destroy(&attr);
    for(t=0; t<threadsNumber; t++) {
       rc = pthread_join(threads[t], NULL);
       if (rc) {
          printf("ERROR; return code from pthread_join() is %d\n", rc);
          exit(-1);
       }
       printf("Main: completed join with thread %ld\n",t);
    }
}

void help(char const *exec, char const opt, char const *optarg) {
	FILE *out = stdout;
	if (opt != 0) {
		out = stderr;
		if (optarg) {
			fprintf(out, "Invalid parameter - %c %s\n", opt, optarg);
		} else {
			fprintf(out, "Invalid parameter - %c\n", opt);
		}
	}
	fprintf(out, "Mandelbrot Set Renderer\n\n");
	fprintf(out, "  -x [0;1]         Center of Re[-1.5;0.5] (default=0.5)\n");
	fprintf(out, "  -y [0;1]         Center of Im[-1;1] (default=0.5)\n");
	fprintf(out, "  -s (0;1]         Inverse scaling factor (default=1)\n");
	fprintf(out, "  -r [pixel]       Image resolution (default=1024)\n");
	fprintf(out, "  -i [iterations]  Iterations or max dwell (default=512)\n");
	fprintf(out, "  -c [colours]     Colour map iterations (default=1)\n");
	fprintf(out, "  -b [block dim]   min block dimension for subdivision (default=16)\n");
	fprintf(out, "  -d [subdivison]  subdivision of blocks (default=4)\n");
	fprintf(out, "  -p [threads]     number of concurrent threads (default=4)\n");
	fprintf(out, "  -m mark Mariani-Silver borders\n");
	fprintf(out, "  -t traditional computation (no Mariani-Silver)\n");
	fprintf(out, "%s [options]  <output-bmp>\n", exec);
}

int main( int argc, char *argv[] )
{
	int ret = 0;
	/* Standard Values */
	char *output = NULL; //output image filepath
	bmpImage *image = NULL; //output image
	double x = 0.5, y = 0.5; // coordinates
	double scale = 1; // scaling factor
	unsigned int colourIterations = 1; //how many times the colour gradient is repeated
	bool quiet = false; //output something or not

	resolution = 1024;
	maxDwell = 512;
	blockDim = 16;
	subdivisions = 4;

	/* Parameter parsing... */
	{
		char c;
		while((c = getopt(argc,argv,"x:y:s:r:o:i:c:b:d:p:mthq"))!=-1) {
			switch(c) {
			case 'x':
				x = clampDouble(atof(optarg),0.0,1.0);
				break;
			case 'y':
				y = clampDouble(atof(optarg),0.0,1.0);
				break;
			case 's':
				scale = clampDouble(atof(optarg),0.0,1.0);
				if (scale == 0) scale = 1;
				break;
			case 'r':
				resolution = atoi(optarg);
				break;
			case 'i':
				maxDwell = atoi(optarg);
				break;
			case 'c':
				colourIterations = atoi(optarg);
				break;
			case 'b':
				blockDim = atoi(optarg);
				blockDim = (blockDim  < 4) ? 4 : blockDim;
				break;
			case 'd':
				subdivisions = atoi(optarg);
				subdivisions = (subdivisions  < 2) ? 2 : subdivisions;
				break;
			case 'm':
				markBorders = true;
				break;
			case 't':
				useMarianiSilver = false;
				break;
			case 'p':
				useThreads = atoi(optarg);
				break;
			case 'q':
				quiet = true;
				break;
			case 'o':
				output = calloc(strlen(optarg) + 1, sizeof(char));
				strncpy(output, optarg, strlen(optarg));
				break;
			case 'h':
				help(argv[0],0, NULL);
				goto exit_graceful;
				break;
			default:
				abort();
			}
		}
	}

	if (output == NULL) {
		fprintf(stderr, "Output argument is not optional!\n");
		goto error_exit;
	}

	/* Initialize the colourMap, which assigns each possible dwell Value a colour */
	/* This might look complex but is just eye candy. */
	/* Could be made much simpler, thus less fancy */

	if (initColourMap(maxDwell / colourIterations) != 0) {
		fprintf(stderr, "Could not initialize colour map!\n");
		goto error_exit;
	};

	/* Putting some magic numbers in place...  */

	double const xmin = -3.5 + (2 * 2 * x);
	double const xmax = -1.5 + (2 * 2 * x);
	double const ymin = -3.0 + (2 * 2 * y);
	double const ymax = -1.0 + (2 * 2 * y);
	double const xlen = fabs(xmin - xmax);
	double const ylen = fabs(ymin - ymax);

	double complex cmax = (xmax - (0.5 * (1 - scale) * xlen)) + (ymax - (0.5 * (1 - scale) * ylen)) * I;
	cmin = (xmin + (0.5 * (1 - scale) * xlen)) + (ymin + (0.5 * (1 - scale) * ylen)) * I;
	dc = cmax - cmin;

	/* Output useful informations... */
	if (!quiet) {
		printf("Center:      [%f,%f]\n",x,y);
		printf("Zoom:        %llu%%\n", (unsigned long long) (1/scale) * 100);
		printf("Iterations:  %u\n", maxDwell);
		printf("Window:      Re[%f,%f], Im[%f,%f]\n", creal(cmin), creal(cmax), cimag(cmin), cimag(cmax));
		printf("Output:      %s\n", output);
	}

	image = newBmpImage(resolution, resolution);
	if (image == NULL) {
		fprintf(stderr, "ERROR: could not allocate bmp image space!\n");
		goto error_exit;
	}

	dwellBuffer = malloc(resolution * resolution * sizeof(dwellType));
	if (dwellBuffer == NULL) {
		fprintf(stderr, "ERROR: could not allocate dwell buffer!\n");
		goto error_exit;
	}
	for (unsigned int i = 0; i < resolution * resolution; i++) {
		dwellBuffer[i] = dwellUncomputed;
	}

	//Threads should start here
	initializeWorkers(useThreads);

	/*if (useMarianiSilver) {
		// Scale the blockSize from res up to a subdividable value
		// Number of possible subdivisions:
		unsigned int const numDiv = ceil(logf((double) resolution/blockDim)/logf((double) subdivisions));
		// Calculate a dividable resolution for the blockSize:
		unsigned int const correctedBlockSize = pow(subdivisions,numDiv) * blockDim;
		// Mariani-Silver subdivision algorithm
		marianiSilver(dwellBuffer, 0, 0, correctedBlockSize);
	} else {
		// Traditional Mandelbrot-Set computation or the 'Escape Time' algorithm
		// computeBlock respects the resolution of the image, so we scale the blocks up to
		// a divideable dimension
		unsigned int block = ceil((double) resolution / useThreads);

		for (unsigned int t = 0; t < useThreads; t++) {
			for (unsigned int x = 0; x < useThreads; x++) {
				escapeTime(dwellBuffer, t * block, x * block, block);
			}
		}
	}*/
	
	//Threads should end here

	// Map dwell buffer to image
	for (unsigned int y = 0; y < resolution; y++) {
		for (unsigned int x = 0; x < resolution; x++) {
			image->rawdata[y * resolution + x] = *getDwellColour(x, y, dwellBuffer[y * resolution + x]);
		}
	}

	// Save the Image
	if(saveBmpImage(image, output)) {
		fprintf(stderr, "ERROR: could not save image to %s\n", output);
		goto error_exit;
	}
	
	goto exit_graceful;
error_exit:
	ret = 1;
exit_graceful:
	if(image)
		freeBmpImage(image);
	if(dwellBuffer)
		free(dwellBuffer);
	freeColourMap();
	
	printf("Main: program completed. Exiting.\n");
	pthread_exit(NULL);
	
	return ret;
}
