//
//  main.cpp
//  Block_Code_15_11
//
//  Created by Stefan Węgrzyn on 20/10/2019.
//  Copyright © 2019 Stefan Węgrzyn. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#define PI 3.141592654

using namespace std;

/*
#define ENCODER_OUTPUT_NUM 2
#define ENCODER_REG_NUM 2
#define ENCODER_STATES_NUM 4
#define ECODER_INPUT_STATES_NUM 2
*/
#define ENCODER_OUTPUT_NUM 2
#define ENCODER_REG_NUM 6
#define ENCODER_STATES_NUM 64
#define ECODER_INPUT_STATES_NUM 2

#define MIN_INPUT_VECTOR_LEN    100000   // 10^(5)
#define MAX_INPUT_VECTOR_LEN    10000000 // 10^(7)

typedef struct
{
    int output[ENCODER_OUTPUT_NUM];
    int final_state[ENCODER_REG_NUM];
} tTransition;
/*
tTransition trussGraph[ENCODER_STATES_NUM][ECODER_INPUT_STATES_NUM] =
                                  { { {0, 0,  0, 0}, {1, 1,  1, 0} },
                                    { {1, 1,  0, 0}, {0, 0,  1, 0} },
                                    { {1, 0,  0, 1}, {0, 1,  1, 1} },
                                    { {0, 1,  0, 1}, {1, 0,  1, 1} } };
*/
tTransition trussGraph[ENCODER_STATES_NUM][ECODER_INPUT_STATES_NUM] =
{
    { {0,0,  0,0,0,0,0,0}, {1,1,  1,0,0,0,0,0} },
    { {1,1,  0,0,0,0,0,0}, {0,0,  1,0,0,0,0,0} },
    { {1,0,  0,0,0,0,0,1}, {0,1,  1,0,0,0,0,1} },
    { {0,1,  0,0,0,0,0,1}, {1,0,  1,0,0,0,0,1} },
    { {0,0,  0,0,0,0,1,0}, {1,1,  1,0,0,0,1,0} },
    { {1,1,  0,0,0,0,1,0}, {0,0,  1,0,0,0,1,0} },
    { {1,0,  0,0,0,0,1,1}, {0,1,  1,0,0,0,1,1} },
    { {0,1,  0,0,0,0,1,1}, {1,0,  1,0,0,0,1,1} },
    { {1,1,  0,0,0,1,0,0}, {0,0,  1,0,0,1,0,0} },
    { {0,0,  0,0,0,1,0,0}, {1,1,  1,0,0,1,0,0} },
    { {0,1,  0,0,0,1,0,1}, {1,0,  1,0,0,1,0,1} },
    { {1,0,  0,0,0,1,0,1}, {0,1,  1,0,0,1,0,1} },
    { {1,1,  0,0,0,1,1,0}, {0,0,  1,0,0,1,1,0} },
    { {0,0,  0,0,0,1,1,0}, {1,1,  1,0,0,1,1,0} },
    { {0,1,  0,0,0,1,1,1}, {1,0,  1,0,0,1,1,1} },
    { {1,0,  0,0,0,1,1,1}, {0,1,  1,0,0,1,1,1} },
    { {1,1,  0,0,1,0,0,0}, {0,0,  1,0,1,0,0,0} },
    { {0,0,  0,0,1,0,0,0}, {1,1,  1,0,1,0,0,0} },
    { {0,1,  0,0,1,0,0,1}, {1,0,  1,0,1,0,0,1} },
    { {1,0,  0,0,1,0,0,1}, {0,1,  1,0,1,0,0,1} },
    { {1,1,  0,0,1,0,1,0}, {0,0,  1,0,1,0,1,0} },
    { {0,0,  0,0,1,0,1,0}, {1,1,  1,0,1,0,1,0} },
    { {0,1,  0,0,1,0,1,1}, {1,0,  1,0,1,0,1,1} },
    { {1,0,  0,0,1,0,1,1}, {0,1,  1,0,1,0,1,1} },
    { {0,0,  0,0,1,1,0,0}, {1,1,  1,0,1,1,0,0} },
    { {1,1,  0,0,1,1,0,0}, {0,0,  1,0,1,1,0,0} },
    { {1,0,  0,0,1,1,0,1}, {0,1,  1,0,1,1,0,1} },
    { {0,1,  0,0,1,1,0,1}, {1,0,  1,0,1,1,0,1} },
    { {0,0,  0,0,1,1,1,0}, {1,1,  1,0,1,1,1,0} },
    { {1,1,  0,0,1,1,1,0}, {0,0,  1,0,1,1,1,0} },
    { {1,0,  0,0,1,1,1,1}, {0,1,  1,0,1,1,1,1} },
    { {0,1,  0,0,1,1,1,1}, {1,0,  1,0,1,1,1,1} },
    { {0,1,  0,1,0,0,0,0}, {1,0,  1,1,0,0,0,0} },
    { {1,0,  0,1,0,0,0,0}, {0,1,  1,1,0,0,0,0} },
    { {1,1,  0,1,0,0,0,1}, {0,0,  1,1,0,0,0,1} },
    { {0,0,  0,1,0,0,0,1}, {1,1,  1,1,0,0,0,1} },
    { {0,1,  0,1,0,0,1,0}, {1,0,  1,1,0,0,1,0} },
    { {1,0,  0,1,0,0,1,0}, {0,1,  1,1,0,0,1,0} },
    { {1,1,  0,1,0,0,1,1}, {0,0,  1,1,0,0,1,1} },
    { {0,0,  0,1,0,0,1,1}, {1,1,  1,1,0,0,1,1} },
    { {1,0,  0,1,0,1,0,0}, {0,1,  1,1,0,1,0,0} },
    { {0,1,  0,1,0,1,0,0}, {1,0,  1,1,0,1,0,0} },
    { {0,0,  0,1,0,1,0,1}, {1,1,  1,1,0,1,0,1} },
    { {1,1,  0,1,0,1,0,1}, {0,0,  1,1,0,1,0,1} },
    { {1,0,  0,1,0,1,1,0}, {0,1,  1,1,0,1,1,0} },
    { {0,1,  0,1,0,1,1,0}, {1,0,  1,1,0,1,1,0} },
    { {0,0,  0,1,0,1,1,1}, {1,1,  1,1,0,1,1,1} },
    { {1,1,  0,1,0,1,1,1}, {0,0,  1,1,0,1,1,1} },
    { {1,0,  0,1,1,0,0,0}, {0,1,  1,1,1,0,0,0} },
    { {0,1,  0,1,1,0,0,0}, {1,0,  1,1,1,0,0,0} },
    { {0,0,  0,1,1,0,0,1}, {1,1,  1,1,1,0,0,1} },
    { {1,1,  0,1,1,0,0,1}, {0,0,  1,1,1,0,0,1} },
    { {1,0,  0,1,1,0,1,0}, {0,1,  1,1,1,0,1,0} },
    { {0,1,  0,1,1,0,1,0}, {1,0,  1,1,1,0,1,0} },
    { {0,0,  0,1,1,0,1,1}, {1,1,  1,1,1,0,1,1} },
    { {1,1,  0,1,1,0,1,1}, {0,0,  1,1,1,0,1,1} },
    { {0,1,  0,1,1,1,0,0}, {1,0,  1,1,1,1,0,0} },
    { {1,0,  0,1,1,1,0,0}, {0,1,  1,1,1,1,0,0} },
    { {1,1,  0,1,1,1,0,1}, {0,0,  1,1,1,1,0,1} },
    { {0,0,  0,1,1,1,0,1}, {1,1,  1,1,1,1,0,1} },
    { {0,1,  0,1,1,1,1,0}, {1,0,  1,1,1,1,1,0} },
    { {1,0,  0,1,1,1,1,0}, {0,1,  1,1,1,1,1,0} },
    { {1,1,  0,1,1,1,1,1}, {0,0,  1,1,1,1,1,1} },
    { {0,0,  0,1,1,1,1,1}, {1,1,  1,1,1,1,1,1} }
};

typedef struct
{
    bool valid;
    int cumulative_distance;
    int previous_state_val;
    int bit_val;
} tGraphNode;

const string inputFilePath = "input.txt";
const string outputFilePath = "output.txt";

//************************************************************************

float gauss(float mean, float sigma);
void kanal(float es_n0, long dl_kan, int *wej, float *wyj);

int getStateVal(int* state);
void getState(int val, int* state);
void performTransition(int* init_state, int input_bit, int* output, int* final_state);
int calcHammingDistance(int* state1, int* state2);
void encodeFrame(int* inputVector, int* encodedInputVector, int frameLen);
void decodeFrame(int* encodedInputVector, int* outputVector, int frameLen);

//************************************************************************

int main(int argc, const char * argv[])
{
    srand((unsigned int)time(NULL));
    
    ifstream file_reader(inputFilePath, ios::in);
    if ( !file_reader.is_open() ) {
        cout << "Could not open file: " << inputFilePath << endl;
        return -1;
    }
    ofstream file_writer(outputFilePath, ios::trunc);
    if ( !file_writer.is_open() ) {
        cout << "Could not open file: " << outputFilePath << endl;
        return -1;
    }

    // Acquire input parameters from the user
    
    long inputLen = 0;
    long encodedInputLen = 0;
    int frameLen = 0;
    long frameNum = 0;
    
    float min_Eb_N0 = -5;  // -5
    float max_Eb_N0 = 10;  // 6-11
    float step_Eb_N0 = 0.25; // 0.5
    
    cout << "Determine number of sent bits in one data frame (proposed values: 100, 500, 1000): ";
    cin >> frameLen;

    cout << "Determine total number of sent bits (min " << MIN_INPUT_VECTOR_LEN;
    cout << " max " << MAX_INPUT_VECTOR_LEN << "): ";
    cin >> inputLen;
    inputLen = inputLen - (inputLen % frameLen);   // make inputLen divisible by frame length
    encodedInputLen = inputLen * ENCODER_OUTPUT_NUM;
    frameNum = inputLen / frameLen;
    
    cout << "Determine MINIMAL value of 'Eb/N0'[dB]: ";
    cin >> min_Eb_N0;
    
    cout << "Determine MAXIMAL value of 'Eb/N0'[dB]: ";
    cin >> max_Eb_N0;
    
    cout << "Determine change step of 'Eb/N0'[dB]: ";
    cin >> step_Eb_N0;
    
    // Read input vector
    char bit[4];
    int* inputVector = new int[inputLen];
    for (long idx = 0; idx < inputLen; idx++)
    {
        file_reader.getline(bit,4);     // "bit + \t\n\r"
        inputVector[idx] = atoi(bit);
    }
    
    cout << "Simulation in progress";
    for (float step = min_Eb_N0; step <= max_Eb_N0; step += step_Eb_N0)
    {
        double BER_no_correction = 0, BER_with_correction = 0;
        
        // NO CORRECTION
     
        // no additional control bits added to the input vector
        
        float* noisyInput = new float[inputLen];
        kanal(step, inputLen, inputVector, noisyInput); // add noise to input bit stream

        for (long idx = 0; idx < inputLen; idx++)
        {
            if ( (noisyInput[idx] >= 0) && (inputVector[idx] == 0) ) BER_no_correction++;
            else if  ( (noisyInput[idx] < 0) && (inputVector[idx] == 1) ) BER_no_correction++;
        }
        
        delete[] noisyInput;
        
        // WITH CORRECTION
        
        // encode every frame separately
        int* encodedInputVector = new int[encodedInputLen];
        for (long frameIdx = 0; frameIdx < frameNum; frameIdx++)
        {
            // determine encoded bit stream for given frame
            encodeFrame(&inputVector[frameIdx*frameLen],
                        &encodedInputVector[frameIdx*frameLen*ENCODER_OUTPUT_NUM],
                        frameLen);
        }
        
        // add noise to encoded input bit stream
        float* noisyEncodedInput_Diff = new float[encodedInputLen];
        kanal(step, encodedInputLen, encodedInputVector, noisyEncodedInput_Diff);
        
        delete[] encodedInputVector;
        
        // convert from differential to single-ended
        int* noisyEncodedInput = new int[encodedInputLen];
        for (long idx = 0; idx < encodedInputLen; idx++)
        {
            if (noisyEncodedInput_Diff[idx] >= 0) noisyEncodedInput[idx] = 1;
            else noisyEncodedInput[idx] = 0;
        }
        
        delete[] noisyEncodedInput_Diff;
        
        // decode every noisy frame separately
        int* outputVector = new int[inputLen];
        for (long frameIdx = 0; frameIdx < frameNum; frameIdx++)
        {
            // determine original bit stream stored in frame
            decodeFrame(&noisyEncodedInput[frameIdx*frameLen*ENCODER_OUTPUT_NUM],
                        &outputVector[frameIdx*frameLen],
                        frameLen);
        }
        
        delete[] noisyEncodedInput;
        
        // determine BER
        for (long idx = 0; idx < inputLen; idx++)
        {
            if ( outputVector[idx] != inputVector[idx] ) BER_with_correction++;
        }
        
        file_writer << "Eb/N0= ";
        file_writer << fixed << setprecision(2) << setw(5) << setfill(' ') << step;
        file_writer << " BER_no_correction= ";
        file_writer << fixed << setprecision(10) << setw(12) << setfill(' ') << BER_no_correction/inputLen;
        file_writer << " BER_with_correction= ";
        file_writer << fixed << setprecision(10) << setw(12) << setfill(' ') << BER_with_correction/inputLen;
        file_writer << endl;
        
        delete[] outputVector;
        
        cout << ".";
    }
    
    cout << endl << "Simulation finished!" << endl;
    delete[] inputVector;
    file_reader.close();
    file_writer.close();
    
    return 0;
}




//******************************************************************

// Function kanal changes binary values into bipolar ones (-1/+1) and adds noise
// *wej - Input vector of binary values (0/1)
// *wyj - Output vector of real numbers
// es_n0 - Es/N0
// dl_kan - the number of input bits
void kanal(float es_n0, long dl_kan, int *wej, float *wyj)
{
    float mean=0;
    float es=1;
    float sygnal;
    float sigma;
    float s_n;
    long y;
    
    s_n=(float) pow(10, (es_n0/10));
    sigma=(float) sqrt (es/(2*s_n));
    
    for (y=0; y<dl_kan; y++)
    {
        sygnal = 2 * *(wej+y)-1; // change the binary value (0/1) into symbol (-1/+1)
        *(wyj+y)=sygnal+gauss(mean,sigma);  // noise addition
    }
}

//*******************************************************************

float gauss(float mean, float sigma)
{
    double x;
    double z;
    
    z=(double)rand()/RAND_MAX;
    if (z==1.0) z=0.9999999;
    x=sigma*sqrt(2.0*log(1.0/(1.0-z)));
    
    z=(double)rand()/RAND_MAX;
    if (z==1.0) z=0.9999999;
    return((float)(mean+x*cos(2*PI*z)));
}

//*******************************************************************

int getStateVal(int* state)
{
    int state_val = 0;
    int powerOfTwo = 1;
    for (int i = 0; i < ENCODER_REG_NUM; i++)
    {
        state_val += state[(ENCODER_REG_NUM-1)-i] * powerOfTwo;
        powerOfTwo *= 2;
    }
    return state_val;
}

//*******************************************************************

void getState(int val, int* state)
{
    int mask = 0b1 << (ENCODER_REG_NUM-1);
    for (int i = 0; i < ENCODER_REG_NUM; i++)
    {
        state[i] = (val & mask) ? 1 : 0;
        mask = mask >> 1;
    }
}

//*******************************************************************

void performTransition(int* init_state, int input_bit, int* output, int* final_state) {
    
    int init_state_val = 0;
    int i = 0;
    
    init_state_val = getStateVal(init_state);
    
    for (i = 0; i < ENCODER_OUTPUT_NUM; i++)
    {
        output[i] = trussGraph[init_state_val][input_bit].output[i];
    }
    
    for (i = 0; i < ENCODER_REG_NUM; i++)
    {
        final_state[i] = trussGraph[init_state_val][input_bit].final_state[i];
    }
}

//*******************************************************************

int calcHammingDistance(int* output1, int* output2)
{
    int distance = 0;
    for (int i = 0; i < ENCODER_OUTPUT_NUM; i++)
    {
        if (output1[i] != output2[i])
        {
            distance++;
        }
    }
    return distance;
}

//*******************************************************************

void encodeFrame(int* inputVector, int* encodedInputVector, int frameLen)
{
    int state[ENCODER_REG_NUM] = {0};   // must start from all-zero state
    int final_state[ENCODER_REG_NUM] = {0};
    for (int i = 0; i < frameLen; i ++)
    {
        performTransition(state, inputVector[i], &encodedInputVector[i*ENCODER_OUTPUT_NUM], final_state);
        for (int j = 0; j < ENCODER_REG_NUM; j++)
        {
            state[j] = final_state[j];
        }
    }
}

//*******************************************************************

void decodeFrame(int* encodedInputVector, int* outputVector, int frameLen)
{
    tGraphNode** trellis = new tGraphNode*[frameLen+1]; // bits in frame + 1
    int sample = 0;
    for (sample = 0; sample < (frameLen+1); sample++)
    {
        trellis[sample] = new tGraphNode[ENCODER_STATES_NUM];
        for (int stateVal = 0; stateVal < ENCODER_STATES_NUM; stateVal++)
        {
            // make sure all nodes are invalid at the beginning
            trellis[sample][stateVal].valid = false;
        }
    }
    
    // make only all-zero state valid at the beginning
    trellis[0][0].valid = true;
    // make initial cumulative distance equal to zero
    trellis[0][0].cumulative_distance = 0;
    
    // fill trellis chart
    for (sample = 0; sample < frameLen; sample++)
    {
        for (int stateVal = 0; stateVal < ENCODER_STATES_NUM; stateVal++)
        {
            
            // ivestigate transitions only from nodes that are already valid
            if (trellis[sample][stateVal].valid == true)
            {
                int output[ENCODER_OUTPUT_NUM] = {0};
                int current_state[ENCODER_REG_NUM] = {0};
                int next_state[ENCODER_REG_NUM] = {0};
                int distance = 0;
                int nextStateVal = 0;
                
                getState(stateVal,current_state);
                
                // assume that bit value is 0/1
                for (int assumed_bit_val = 0; assumed_bit_val < ECODER_INPUT_STATES_NUM; assumed_bit_val++)
                {
                    performTransition(current_state, assumed_bit_val, output, next_state);
                    distance = calcHammingDistance(&encodedInputVector[sample*ENCODER_OUTPUT_NUM], output);
                    nextStateVal = getStateVal(next_state);
                    
                    // modify nodes from next sample
                    if (trellis[sample+1][nextStateVal].valid == false)
                    {
                        trellis[sample+1][nextStateVal].valid = true;
                        trellis[sample+1][nextStateVal].cumulative_distance = trellis[sample][stateVal].cumulative_distance + distance;
                        trellis[sample+1][nextStateVal].previous_state_val = stateVal;
                        trellis[sample+1][nextStateVal].bit_val = assumed_bit_val;
                    }
                    else
                    {
                        // if accessing already valid node, update it only if
                        // this transition assures smaller cumulative distance
                        if (trellis[sample][stateVal].cumulative_distance + distance < trellis[sample+1][nextStateVal].cumulative_distance)
                        {
                            trellis[sample+1][nextStateVal].cumulative_distance = trellis[sample][stateVal].cumulative_distance + distance;
                            trellis[sample+1][nextStateVal].previous_state_val = stateVal;
                            trellis[sample+1][nextStateVal].bit_val = assumed_bit_val;
                        }
                    }
                }
            }
        }
    }
    
    // find most probable path
    int most_probable_final_state_val = 0;
    int smallest_distance = trellis[frameLen][0].cumulative_distance;
    for (int stateVal = 1; stateVal < ENCODER_STATES_NUM; stateVal++)
    {
        // search in final sample for smallest cumulative distance
        if (trellis[frameLen][stateVal].valid == true)
        {
            if (trellis[frameLen][stateVal].cumulative_distance < smallest_distance)
            {
                smallest_distance = trellis[frameLen][stateVal].cumulative_distance;
                most_probable_final_state_val = stateVal;
            }
        }
    }
    
    // determine most probable frame content
    int chosen_state = most_probable_final_state_val;
    for (sample = frameLen; sample > 0; sample--)
    {
        outputVector[sample-1] = trellis[sample][chosen_state].bit_val;
        chosen_state = trellis[sample][chosen_state].previous_state_val;
    }
    
    for (sample = 0; sample < frameLen; sample++)
    {
        delete[] trellis[sample];
    }
    delete[] trellis;
}

//*******************************************************************
