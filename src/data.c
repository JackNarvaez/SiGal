#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 

#define SETUP_NAME_SIZE 32  

void read_parameters(const char *File_address, double prmts[], char *setup)
{
    FILE *File;
    File = fopen(File_address, "r"); 
    char line[256];
    int row = 0;
    while (fgets(line, sizeof(line), File)) {  
        if (line[0] == '\n' || line[0] == '#')  
            continue;
        
        char *token = strtok(line, "\t ");  
        if (row < 6 && token != NULL) {  
            prmts[row] = atof(token);  
            row += 1;
        } else if (token != NULL) { 
            // Copy the setup name into the setup variable
            strncpy(setup, token, SETUP_NAME_SIZE - 1);
            setup[SETUP_NAME_SIZE - 1] = '\0'; 
            for (int i = strlen(setup) - 1; i >= 0 && isspace(setup[i]); i--) {
                setup[i] = '\0';
            }
        }
    }
    fclose(File); 
}

void read_data(const char *File_address, double *Pos, double *Vel, double *Mass) {
    /*---------------------------------------------------------------------------
    Reads data from bodies: position, velocity and mass.
    -----------------------------------------------------------------------------
    Arguments:
      File_address:   File address from where the data is read.
      Pos         :   Position (1D vector) [xi, yi, zi].
      Vel         :   Velocity (1D vector) [vxi, vyi, vzi].
      Mass        :   Mass (1D vector) [mi].
    ---------------------------------------------------------------------------*/

    FILE *File;
    File = fopen(File_address, "r");
    char line[256];
    int row = 0;
    int ii;
    while (fgets(line, sizeof(line), File)) {
        if (line[0] == '\n' || line[0] == '#')
            continue;
        else {
            char *token;
            token = strtok(line, "\t");
            for (ii = 0; ii < 3; ii++) {
                Pos[ii + 3 * row] = atof(token);
                token = strtok(NULL, "\t");
            }
            for (ii = 0; ii < 3; ii++) {
                Vel[ii + 3 * row] = atof(token);
                token = strtok(NULL, "\t");
            }
            Mass[row] = atof(token);
            row += 1;
        }
    }
    fclose(File);
}

void write_data(const char *File_address, double *Pos, double *Vel, double *Mass, int N) {
    /*---------------------------------------------------------------------------
    Writes data about bodies: position, velocity and mass.
    -----------------------------------------------------------------------------
    Arguments:
      File_address:   File address from where the data is read.
      Pos         :   Position (1D vector) [xi, yi, zi].
      Vel         :   Velocity (1D vector) [vxi, vyi, vzi].
      Mass        :   Mass (1D vector) [mi].
      N           :   Number of bodies
    ---------------------------------------------------------------------------*/

    FILE *File;
    File = fopen(File_address, "w");
    int ii;
    fprintf(File, "#\t%d\n", N);
    for (ii=0; ii<N; ii++) {
        fprintf(File, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Pos[3*ii], Pos[3*ii+1], Pos[3*ii+2], Vel[3*ii], Vel[3*ii+1], Vel[3*ii+2], Mass[ii]);
    }
    fclose(File);
}