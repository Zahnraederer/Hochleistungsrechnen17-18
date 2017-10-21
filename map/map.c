#include <stdio.h>

// Definieren Sie ein enum cardd

typedef enum {N=1, E=2, S=4, W=8} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält

cardd map[3][3]; /*= { {0,N,0},
                    {W,0,E},
                    {0,S,0}};
*/

char give_char(int pos, char* peter)
{
  char x;
  switch(pos) {
    case 0   : x = '0'; break;
    case 1   : x = 'N'; break;
    case 2   : x = 'E'; break;
    case 3   : *peter = 'E'; x = 'N'; break;
    case 4   : x = 'S'; break;
    case 8   : x = 'W'; break;
    default  : x = 'X';
  }
  return x;
}

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  //if(dir == 3)
    //printf("%d  %d  %d\n",x,y,dir);
  map[x][y] = dir;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  char peter = ' ';
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      printf("%c%c  ",give_char(map[i][j],&peter),peter);
      peter = ' ';
      //printf("%d   ",map[i][j]);
    }
    printf("\n");
  }
  //printf("%d",N|E);
}


int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	//set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	//set_dir(2, 2, E|W);

	show_map();

	return 0;
}
