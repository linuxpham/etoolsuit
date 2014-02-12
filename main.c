/*******************************************************************************
 * File:   main.c
 * Author: Linuxpham
 * Description : Etool-Suit for Biologist
 * Created on February 11, 2014, 11:00 PM
 ******************************************************************************/
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <stddef.h>

#ifdef _WIN64
#include <windows.h>
#define SPLIT_PATH_NAME "\\"
#define SPLIT_PATH_NAME_CHAR '\\'
#define _TYPE_CMD 1
#define _IS_WINDOW 1
#define unlink _unlink
#define mkdir _mkdir
#elif _WIN32
#include <windows.h>
#define SPLIT_PATH_NAME "\\"
#define SPLIT_PATH_NAME_CHAR '\\'
#define _TYPE_CMD 1
#define _IS_WINDOW 1
#define unlink _unlink
#define mkdir _mkdir
#else
#define SPLIT_PATH_NAME "/"
#define SPLIT_PATH_NAME_CHAR '/'
#ifdef __APPLE__
#define _TYPE_CMD 2
#define _OS_PLATFORM_ID 1
#else
#define _TYPE_CMD 3
#define _OS_PLATFORM_ID 2
#endif
#endif
#define _VERSION "1.0.10"
#define _SUIT_NAME "suit"
#define _BWA_NAME "bwa"
#define _SAM_NAME "sam"
#define _BED_NAME "bed"
#define _IGV_NAME "igv"
#define _BAM_NAME "bam"
#define _SIBELIA_NAME "sib"
#define _SYS_INI "sys.ini"
#define _MAX_BUFFER 255
#define _OUT_PREPARE_NAME "out"

/**
 * Create directory
 * @param sDirectoryPath
 * @return integer
 */
int isMkDirectory(char *sDirectoryPath) {
    //Check access directory status
    if (access(sDirectoryPath, 0) == 0) {
        //Get directory status information
        struct stat arrStatus;
        stat(sDirectoryPath, &arrStatus);

        //Check directory exist
        if (arrStatus.st_mode & S_IFDIR) {
            return 1;
        }
        return -1;
    } else {
        int iResponse = -1;
#if defined(_WIN32)
        iResponse = mkdir(sDirectoryPath);
#else
        iResponse = mkdir(sDirectoryPath, 0755);
#endif
        if (iResponse == 0) {
            return 1;
        }
        return 0;
    }
    return 0;
}

/**
 * Show help table
 */
static void showHelp(void) {
    fprintf(stdout, "\nNAME"
            "\n\tEtool Suit Application for Biologist"
            "\n\nVERSION"
            "\n\t%s"
            "\n\nSYNOPSIS"
            "\n\tetool\t<command>"
            "\n\t\t+ pre <reference fasta> <left fastQ genome> <right fastQ genome> <output directory> <thread number>: Prepare processing"
            "\n\t\t+ bwa <command> [options]: BWA tool"
            "\n\t\t+ sam <command> [options]: Sam tool"
            "\n\t\t+ bed <command> [options]: Bed tool"
            "\n\t\t+ bmv : Bam viewer"
            "\n\t\t+ igv : IGV browser"
            "\n\t\t+ sib <options>: Sibelia"
            "\n\t\t+ csb <options>: C-Sibelia", _VERSION);
    fprintf(stdout, "\n\nEXAMPLES"
            "\n\tetool pre reference.fasta left.fastq.gz right.fastq.gz /tmp/test 4 (Prepare processing)"
            "\n\tetool bwa index reference.fasta (Index reference genome)"
            "\n\tetool bwa mem -a -t 10 reference.fasta left.fastq.gz right.fastq.gz > out.sam (FastQ alignment)"
            "\n\tetool sam view -bS out.sam > out.bam (Create BAM)"
            "\n\tetool sam sort out.bam out.sorted (Sort BAM)"
            "\n\tetool sam index out.sorted.bam (Index sorted BAM)"
            "\n\tetool bed bamtobed -i out.sorted.bam > out.bed (BAM to BED)"
            "\n\tetool bmv (BAM viewer)"
            "\n\tetool igv (IGV browser)"
            "\n\tetool sib (Sibelia)"
            "\n\tetool csb (C-Sibelia)");
    fprintf(stdout, "\n\nDIAGNOSTICS"
            "\n\tThe etool utility exits 0 on success, and >0 if an error occurs.");
    fprintf(stdout, "\n");
}

/**
 * Check exist command
 * @param sCommand
 * @return integer
 */
static int isExist(char *sCommand) {
    return (access(sCommand, F_OK) == 0);
}

/**
 * Main application
 * @param argc
 * @param argv
 * @return <int>
 */
int main(int argc, char** argv) {
    if (argc < 2) {
        showHelp();
        return (EXIT_FAILURE);
    }

    //Create command and parameters
    char sCommandProc[30] = {'\0'};
    char sDirectorySub[30] = {'\0'};

    //If use BAM tool
    if (strcmp(argv[1], "bwa") == 0) {
        strcpy(sCommandProc, "bwa");
        strcpy(sDirectorySub, "bwa");
    }//If use SAM tool
    else if (strcmp(argv[1], "sam") == 0) {
        strcpy(sCommandProc, "samtools");
        strcpy(sDirectorySub, "sam");
    }//If use BED tool
    else if (strcmp(argv[1], "bed") == 0) {
        strcpy(sCommandProc, "bedtools");
        strcpy(sDirectorySub, "bed");
        strcat(sDirectorySub, SPLIT_PATH_NAME);
        strcat(sDirectorySub, "bin");
    } else if (strcmp(argv[1], "bmv") == 0) {
        switch (_TYPE_CMD) {
            case 1:
                strcpy(sCommandProc, "BamView.bat");
                break;
            case 3:
                strcpy(sCommandProc, "BamView.sh");
                break;
            case 2:
                strcpy(sCommandProc, "BamView.command");
                break;
            default:
                break;
        }
        strcpy(sDirectorySub, "bam");
    }//If use IGV tool
    else if (strcmp(argv[1], "igv") == 0) {
        switch (_TYPE_CMD) {
            case 1:
                strcpy(sCommandProc, "igv.bat");
                break;
            case 3:
                strcpy(sCommandProc, "igv.sh");
                break;
            case 2:
                strcpy(sCommandProc, "igv.command");
                break;
            default:
                break;
        }
        strcpy(sDirectorySub, "igv");
    }//If use SIBELIA tool
    else if (strcmp(argv[1], "sib") == 0) {
        strcpy(sCommandProc, "Sibelia");
        strcpy(sDirectorySub, "sib");
        strcat(sDirectorySub, SPLIT_PATH_NAME);
        strcat(sDirectorySub, "bin");
    }//If use CSIBELIA tool
    else if (strcmp(argv[1], "csb") == 0) {
        strcpy(sCommandProc, "C-Sibelia.py");
        strcpy(sDirectorySub, "sib");
        strcat(sDirectorySub, SPLIT_PATH_NAME);
        strcat(sDirectorySub, "bin");
    } else if (strcmp(argv[1], "pre") == 0) {
        strcpy(sCommandProc, "pre");
        strcpy(sDirectorySub, "pre");
    }

    //Check command length data
    if (strlen(sCommandProc) == 0) {
        showHelp();
        return (EXIT_FAILURE);
    }

    //Get current directory path
    char sDirectoryProc[strlen(argv[0]) + 1];
    memset(&sDirectoryProc, '\0', sizeof (sDirectoryProc));
    strcpy(sDirectoryProc, argv[0]);
    char *cPosition = strrchr(sDirectoryProc, SPLIT_PATH_NAME_CHAR);
    if (cPosition == NULL) {
        showHelp();
        return (EXIT_FAILURE);
    }
    memset(&sDirectoryProc, '\0', sizeof (sDirectoryProc));
    strncpy(sDirectoryProc, argv[0], (cPosition - sDirectoryProc + 1));

#ifdef DEBUG
    fprintf(stdout, "Current working directory: %s\n", sDirectoryProc);
#endif

    //Create command parameters
    char *sCommandParameters;
    int iLoop = 0;

    //Release command parameters
    sCommandParameters = (char *) calloc(1, sizeof (char));
    if (sCommandParameters == NULL) {
        showHelp();
        return (EXIT_FAILURE);
    }

    //Loop to put parameters data
    for (iLoop = 2; iLoop < argc; iLoop++) {
        sCommandParameters = realloc(sCommandParameters, strlen(sCommandParameters) + 10 + strlen(argv[iLoop]));
        strcat(sCommandParameters, argv[iLoop]);
        if (iLoop < (argc - 1)) {
            strcat(sCommandParameters, " ");
        }
    }

    //Create full command to execute
    char sCommandFull[strlen(sCommandParameters) + 512];
    memset(sCommandFull, '\0', sizeof (sCommandFull));

    //Check full command to execute
    if (strcmp(argv[1], "pre") == 0) {
        if (argc < 5) {
            showHelp();
            return (EXIT_FAILURE);
        }
                
        //Create command
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bwa");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bwa");

        //Check command is existing
        if (!isExist(sCommandFull)) {
            if (sCommandParameters) {
                free(sCommandParameters);
                sCommandParameters = NULL;
            }
            showHelp();
            return (EXIT_FAILURE);
        }

        //Index reference file
        strcat(sCommandFull, " index ");
        strcat(sCommandFull, argv[2]);

        //Create threads number
        char outputDir[_MAX_BUFFER] = {'\0'};
        if (argv[5] == NULL) {
            strcpy(outputDir, "./");
        } else {
            strncpy(outputDir, argv[5], (_MAX_BUFFER - 1));

            if (isMkDirectory(outputDir) == 0) {
                showHelp();
                return (EXIT_FAILURE);
            }
        }

        //Create threads number
        char threadNum[10] = {'\0'};
        if (argv[6] == NULL) {
            strcpy(threadNum, "4");
        } else {
            strncpy(threadNum, argv[6], 9);
        }

        //Build SAM file
        strcat(sCommandFull, " && ");
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bwa");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bwa mem -a -t ");
        strcat(sCommandFull, threadNum);
        strcat(sCommandFull, " ");
        strcat(sCommandFull, argv[2]);
        strcat(sCommandFull, " ");
        strcat(sCommandFull, argv[3]);
        strcat(sCommandFull, " ");
        strcat(sCommandFull, argv[4]);
        strcat(sCommandFull, " > ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".sam");

        //Convert from SAM to BAM
        strcat(sCommandFull, " && ");
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "sam");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "samtools view -bS ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".sam");
        strcat(sCommandFull, " > ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".bam");

        //Sort BAM file
        strcat(sCommandFull, " && ");
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "sam");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "samtools sort ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".bam");
        strcat(sCommandFull, " ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".sorted");

        //Index SORTED BAM file
        strcat(sCommandFull, " && ");
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "sam");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "samtools index ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".sorted.bam");

        //Convert from BAM to BED file
        strcat(sCommandFull, " && ");
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bed");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bin");
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, "bedtools bamtobed -i ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".sorted.bam");
        strcat(sCommandFull, " > ");
        strcat(sCommandFull, outputDir);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, _OUT_PREPARE_NAME);
        strcat(sCommandFull, ".bed");
    } else {
        strcat(sCommandFull, sDirectoryProc);
        strcat(sCommandFull, _SUIT_NAME);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, sDirectorySub);
        strcat(sCommandFull, SPLIT_PATH_NAME);
        strcat(sCommandFull, sCommandProc);

        //Check command is existing
        if (!isExist(sCommandFull)) {
            if (sCommandParameters) {
                free(sCommandParameters);
                sCommandParameters = NULL;
            }
            showHelp();
            return (EXIT_FAILURE);
        }

        //Add extended command data
        strcat(sCommandFull, " ");
        strcat(sCommandFull, sCommandParameters);
    }

#ifdef DEBUG
    fprintf(stdout, "Use the [%s] command\n", sCommandFull);
#endif    

    //Execute the full command    
    FILE *outputStream = popen(sCommandFull, "r");
    char arrBuffer[_MAX_BUFFER] = {'\0'};
    while (fgets(arrBuffer, _MAX_BUFFER, outputStream) != NULL) {
        fprintf(stdout, "%s", arrBuffer);
        memset(&arrBuffer, '\0', sizeof (arrBuffer));
    }
    pclose(outputStream);

    //Free memory leaking
    if (sCommandParameters) {
        free(sCommandParameters);
        sCommandParameters = NULL;
    }

    //Return default data
    return (EXIT_SUCCESS);
}

