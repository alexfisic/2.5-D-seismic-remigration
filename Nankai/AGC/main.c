/******************************************************************************
**
** Function:         : main()
**
** Description       : This files reads the 16bit Signed data from a file 
**                     convert it into -1 to +1 range and calls the AGC 
**                     AGC function and again the output data in the range
**                     of -1 to +1 converts in to 16bit signed data and 
**                     write back in to the file.
**
** Programmer        : Jaydeep Appasaheb Dhole
**                   : Associate Software Engineer ( DSP )
**                   : Aparoksha Tech. Pvt. Ltd. ,Bangalore,INDIA.
**                   : http://www.aparoksha.com
**                   : <jaydeepdhole@gmail.com>
**
**
** Date              : 26 May 2006.
**
******************************************************************************/

#include<stdio.h>
#include<math.h>

#define L_FRAME      256         /* number of samples or frame length */
#define GAIN_LEVEL   -30         /* output power in db                */

void AGC(float *, float *, float, int);

void main(void)
{
    int      i,k;
    float    input_buff[L_FRAME];
    float    output_buff[L_FRAME];
    int      output_buff_16bit[L_FRAME];
    float    energy;
    double   P;
    FILE     *fs,*fd;

    fs = fopen("input.dat","r");
    if(fs==NULL)
        printf("\n Error in opening input file...");

    fd = fopen("output.dat","w");
    if(fd==NULL)
        printf("\n Error in opening output file...");

    /* Read the 16 bit input data from the file */
    /* convert it in to the range of -1 to +1 */
    for(i = 0; i < L_FRAME; i++)
    {
        fscanf(fs,"%d",&k);
        input_buff[i] = (float)k/32767;
    }

    /* input energy Calculation */
    energy = 0;
    for(i = 0; i < L_FRAME; i++)
        energy += input_buff[i] * input_buff[i];
    P = 10 * log10(energy/L_FRAME);
        printf("\n input power is %f\n",(float)P);

    AGC(input_buff,output_buff,GAIN_LEVEL,L_FRAME);

    /* output energy Calculation */
    energy = 0;
    for(i = 0; i < L_FRAME; i++)
        energy += output_buff[i] * output_buff[i];
    P = 10 * log10(energy/L_FRAME);
        printf("\n output power in floating point arithmetic is %f\n",(float)P);

    /* Output in the range of -1 to +1                  */
    /* Convert it in to -32768 to +32767 i.e.16 bit format */
    /* Make it round by adding 0.5 to reduce the error  */
    for(i = 0; i < L_FRAME; i++)
        output_buff_16bit[i] = (int)(output_buff[i] * 32767 + 0.5);

    /* Save the 16 bit output into file */
    for(i = 0; i < L_FRAME; i++)
        fprintf(fd,"%d\n",output_buff_16bit[i]);

    fclose(fs);
    fclose(fd);
}