#include <stdio.h>		// Bibliothèque standard d'entrées/sorties (printf, fprintf, fopen, etc.)
#include <stdlib.h>		// Bibliothèque standard (exit, malloc, free, etc. si besoin)
#include <math.h>		// Pour les fonctions mathématiques (sqrt, cos, sin, fabs, etc.)


/////////////////////////////////////////////////////// Équation différentielle des géodésiques ///////////////////////////////////////////////////////
// u'' = 3/2 u^2 - u 
double du_point_dtheta(double u) {
	return ((3.0/2.0) * u * u) - u; // Notre équation de trajectoire
}



///////////////////////////////////////////////////////// Schéma de Runge-Kutta d'ordre 4 /////////////////////////////////////////////////////////////
// Runge Kutta4
void RK4(FILE *file, double u, double u_point, double theta, double dtheta, int Ntheta, double rs, int color, int resolution) {
    
	double r=rs/u, x=r*cos(theta), y=r*sin(theta); // On calcule ce qui nous intéresse pour le futur plot

    fprintf(file, "%.15e \t\t %.15e \t %d\n", x, y, color); // On remplit la 1ere ligne avec les 1eres valeurs

	for (int i=0; i<Ntheta; i++) {
		double k1 = (u_point) * dtheta;								// Coefficients RK4 k1  pour u
		double k1_point = du_point_dtheta(u) * dtheta;				// Coefficients RK4 k1' pour u'

		double k2 = (u_point + k1_point/2.0) * dtheta;				// Coefficients RK4 k2  pour u
		double k2_point = du_point_dtheta(u + k1/2.0) * dtheta;		// Coefficients RK4 k2' pour u'

		double k3 = (u_point + k2_point/2.0) * dtheta;				// Coefficients RK4 k3  pour u
		double k3_point = du_point_dtheta(u + k2/2.0) * dtheta;		// Coefficients RK4 k3' pour u'

		double k4 = (u_point + k3_point) * dtheta;					// Coefficients RK4 k4  pour u
		double k4_point = du_point_dtheta(u + k3) * dtheta;			// Coefficients RK4 k4' pour u'
		
		u = u + (k1 + 2*k2 + 2*k3 + k4)/6.0; 										// Nouvelle valeur de u
		u_point = u_point + (k1_point + 2*k2_point + 2*k3_point + k4_point)/6.0; 	// Nouvelle valeur de u'
        theta += dtheta; 															// On avance l’angle θ de dθ

        r=rs/u;
        x=r*cos(theta);
        y=r*sin(theta);
        if (r<rs || r>1e8) break; 		// Si le rayon tombe dans le trou noir ou qu'il est beaucoup trop loin, alors on arrête le calcul
        if ((i%resolution)==0) { 		// Mettre le modulo n (i%n) qu'on souhaite pour diminuer la taille du fichier (la résolution finale de l'image) sans changer la précision
            fprintf(file, "%.15e \t\t %.15e \t %d\n", x, y, color); //On remplit ce qu'on vient de calculer
        }
	}
}


static void progress_bar(int current, int total) {
    const int width = 50; 							 	// largeur de la barre
    double ratio = (double)current / (double)total;		// Fraction du travail accompli (entre 0 et 1)
    int filled = (int)(ratio * width); 					// Nombre de caractères "remplis" dans la barre

    printf("\r[");										// '\r' ramène le curseur au début de la ligne pour réécrire dessus
    for (int i = 0; i < width; ++i)
        putchar(i < filled ? '=' : ' ');				// On met '=' pour la partie faite, espace pour le reste
    printf("] %3d%%", (int)(ratio * 100.0));			// Affiche le pourcentage sous forme " xx%"
    fflush(stdout);  									// => La barre se remplace elle même au final
}


void explorationCI(double uc_max, double uc_min, double theta_min, double u0, FILE *file, double dtheta, int Ntheta, double rs, int resolution) {
	fprintf(file, "\t\t x \t\t\t\t\t\t y \t\t\t\t\t\t color\n"); 		// En-tête de mon .txt
	for (double uc=uc_max, step=0; uc>=uc_min; uc-=uc_min, step++) {
	
		double theta = theta_min;												// On initialise l'angle θ au début de la trajectoire
        double u = u0; 															// On fixe la valeur initiale de u (rayon très grand ⇒ u0 très petit)
        double u_point = sqrt( (4.0/27.0) * uc*uc - u0*u0 + u0*u0*u0);			// Condition initiale sur u'
		
		int color; double epsilon = 1e-6;										// On définit une "couleur" (un code entier) en fonction de la valeur de uc
		if (uc > 1.0 + epsilon) color = 0;										// Rayon "non critique" (au-dessus de la valeur critique)
		else if (fabs(uc - 1.0) < epsilon) color = 1;							// Rayon "critique" (uc ~ 1)
		else color = 2;															// Rayon "en dessous du critique"
		
		RK4(file, u, u_point, theta, dtheta, Ntheta, rs, color, resolution);	// On appelle notre intégrateur RK4 pour tracer la géodésique

		fprintf(file, "\t\t NaN \t\t\t\t\t\t NaN \t NaN\n");					// Permet de dire à python qu'on va changer de rayon
		
		
		progress_bar(step, (int)((uc_max - uc_min) / uc_min));					// Mise à jour de la barre de progression
	}
    printf("\n");		// On passe à la ligne après la barre de progression
	fclose(file);		// On ferme le fichier une fois fini
}




int main(){
	// Ouverture du fichier pour stocker les valeurse
	FILE* file = fopen("eq_diff_1.txt", "w");

	double G = 6.67430e-11;  					// Constante gravitationnelle
	double c = 2.99792458e8; 					// Vitesse de la lumière
	double M_sun = 1.98847e30; 					//Masse du soleil
	double masse_trou_noir = 5 * M_sun; 		// Masse du trou noir
	double rs = 2 * G * masse_trou_noir/(c*c);  // Rayon de Schwartzchilds
	

	double u0 = 1e-3;  							// r0 = rs/u0;
	double theta_max = 100.0, theta_min = 0;	// valeur de theta initiale et finale
	double dtheta=1e-6;							// Pas de theta
	
	double uc_min = 0.01; 						//Valeur finale de uc (tend vers 0)
	double uc_max = 100;   						//Valeur initiale de uc (tend vers l'infini)
	int resolution = 5000;						//Resolution du au modulo dans RK4
	
	int Ntheta = (theta_max-theta_min)/dtheta;	// Nombre d'élément à calculer
    
    explorationCI(uc_max, uc_min, theta_min, u0, file, dtheta, Ntheta, rs, resolution);	// On lance l'exploration de toutes les conditions initiales en uc
    		
	return 0; 	// Fin du programme principal : tout s'est bien passé, on retourne 0.
}

