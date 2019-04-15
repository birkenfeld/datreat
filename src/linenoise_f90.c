/*
   Command Parser interface to linenoise library

   Author: Piotr Zolnierczuk (FZ  Juelich)
   Date:


 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "linenoise.h"

#define MAX_COMMAND 256

static char* historyFileName=NULL;
static char* completionList[MAX_COMMAND] = { NULL, };
static int   completionListSize = 0;

/*
   The completion function
 */
static void
fclinenoise_complete(const char *prefix, linenoiseCompletions *lc) {
   int n = strlen(prefix);
   int k = 0;
   while(k<completionListSize && completionList[k]!=NULL) {
     if ( !strncmp(prefix, completionList[k], n))
       linenoiseAddCompletion(lc,completionList[k]);
     k++;
   }
}


/*
   Initialize linenoise library, set up completion and history
 */
int
fclinenoise_init(char *filename)
{
  int filenameLen=strlen(filename);
  /* Set the completion callback. This will be called every time the
   * user uses the <tab> key. */
  linenoiseSetCompletionCallback(fclinenoise_complete);

  /* Load history from file. The history file is just a plain text file
   * where entries are separated by newlines. */
  linenoiseHistoryLoad(filename); /* Load the history at startup */

  /* save history filename */
  historyFileName = (char *)malloc(filenameLen+1);
  strncpy(historyFileName, filename, filenameLen);
  historyFileName[filenameLen]=0x00;
  return 0;
}


/*
   Add command to the completion list
 */
int
fclinenoise_add(char *command)
{
  int cmdlen = strlen(command);
  char *buffer;
  if ( completionListSize>=MAX_COMMAND )
    return -1;
  buffer   = malloc(cmdlen+1);
  strncpy(buffer, command, cmdlen+1);
  completionList[completionListSize] = buffer;
  completionListSize++;
  return 0;
}


/*
   Read line
 */
int
fclinenoise_read(int mylen, char *myline, char* prompt)
{
  char *line;
  int   len;
  int   res=0;

  line = linenoise(prompt);
  if (line==NULL)
    return -1;
  len = 0;
  memset(myline, 0x20, mylen);   /* Set it to all spaces */
  if (line[0] != '\0' && line[0] != '/') {
    len=strlen(line);
    linenoiseHistoryAdd(line);             /* add to the history. */
    linenoiseHistorySave(historyFileName); /* save the history on disk. */
    strncpy(myline, line, mylen<len ? mylen : len);
  } else if (!strncmp(line,"/historylen",11)) {
    len = atoi(line+11);
    linenoiseHistorySetMaxLen(len);
  } else if (!strncmp(line,"/clear",6)) {
    linenoiseClearScreen();
  }
  //exit:
  free(line);
  return res;
}
