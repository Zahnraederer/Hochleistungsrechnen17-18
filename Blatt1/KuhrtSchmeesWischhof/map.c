#include <stdio.h>

// Definieren Sie ein enum cardd

typedef enum {N=1, E=2, S=4, W=8} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält

cardd map[3][3]; /*= { {0,N,0},
                    {W,0,E},
                    {0,S,0}};
*/

//Hilfsfunktion die mithilfe von Switch den Binärwerten die entsprechende
//Zeichenkette zuweist
char* give_char(int pos)
{
  char* x;
  switch(pos) {
    case 0   : x = "0"; break;
    case N   : x = "N"; break;
    case E   : x = "E"; break;
    case S   : x = "S"; break;
    case W   : x = "W"; break;

    case N|E   : x = "NE"; break;
    case S|E   : x = "SE"; break;
    case S|W   : x = "SW"; break;
    case N|W   : x = "NW"; break;

    default  : x = "X";
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

  if (x < 0 || x > 2) {
    fprintf(stderr, "x out of bounds: %d when 0-2 is expected\n", x);
  } else if (y < 0 || y > 2) {
    fprintf(stderr, "y out of bounds: %d when 0-2 is expected\n", y);
  } else if (dir < 0 /* Untere Begrenzung */
          || dir > (S|W) /* Obere Begrenzung */
          || (dir & (N|S)) == (N|S) /* N- und S-Bit gleichzeitig gesetzt? */
          || (dir & (E|W)) == (E|W)) /* W- und E-Bit gleichzeitig gesetzt? */ {
    fprintf(stderr, "Unknown value for dir: %d\n", dir);
  } else {
    map[y][x] = dir;
  }
}


// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  char* line[3];
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      line[j] = give_char(map[i][j]);
    }
    printf("%-4s%s%4s\n",line[0],line[1],line[2]);
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
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);

	show_map();

	return 0;
}
