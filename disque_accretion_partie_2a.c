#include <stdio.h>
#include <stdlib.h> 
#include <math.h>


////////////////////////////////////////////// Fonction qui permet d'afficher une barre de progression //////////////////////////////////////////////
void progress_bar(double progress){
    const int barWidth = 50;                // Longueur de la barre (nombre de caractères)
    
    if (progress < 0.0) progress = 0.0;     // On s'assure que progress reste dans [0, 1]
    if (progress > 1.0) progress = 1.0;

    int pos = (int)(progress * barWidth);   // Position actuelle dans la barre (index de la "tête" de la progression)

    printf("\r[");                          // '\r' ramène le curseur au début de la même ligne pour réécrire par-dessus
 

    for (int i = 0; i < barWidth; i++) {    // On parcourt chaque "case" de la barre
        if (i < pos)        printf("=");    // Partie déjà complétée : on affiche des '='
        else if (i == pos)  printf(">");    // Position actuelle : on met un '>' pour symboliser l'avancement
        else                printf(" ");    // Partie restante : on laisse des espaces vides
    }

    printf("] %6.2f%%", progress * 100.0);  // Affiche le pourcentage avec 2 décimales, sur 6 caractères de large
    fflush(stdout);                         // On force l'affichage immédiat (évite que ça reste dans le buffer)
}




//////////////////// Fonction qui permet de calculer theta_d en fonction de alpha et i, c'est l'équation qui est dans le sujet ////////////////////
double theta_d(double alpha,double i){
    if(cos(alpha)*cos(alpha)*cos(i)*cos(i) == 1.0){     // On vérifie simplement si le dénominateur est égal à 0
        return acos((-sin(alpha)*cos(i))/sqrt(1e-6));
    }
    return acos((-sin(alpha)*cos(i))/sqrt(1-cos(alpha)*cos(alpha)*cos(i)*cos(i)));
}


void explorationCI(FILE* file, double i, double b_start, double b_end, double db, double alpha_start, double dalpha, double R_min, double R_max){
    fprintf(file, "\t\t\t %-30s %-35s %-40s\n", "b", "alpha", "r");         // Entête de notre fichier .txt (%-30.. pour essayer de centrer sur les colonnes)

    int n_steps = (int)((b_end - b_start) / db);                            // Pour la barre de progression, calcule le nombre total 
    int step = 0;                                                           // On initialise step qui augmentera dans la boucle

    for(double b=b_start; b<b_end ; b+=db, step++)                          // On boucle sur tous les éléments b pour remplir le détecteur
    {
        progress_bar((step + 1) / (double)n_steps);                         // Barre de progression qui se met à jour

        for(double alpha = alpha_start ; alpha<2*M_PI ; alpha+=dalpha)      // On boucle sur tous les angles alpha pour remplir le détecteur sur 2PI
        {   
            double theta_end = theta_d(alpha,i);                            // lorsqu'on touche le disque d'accrétion, c'est qu'on a atteint l'angle theta_d, donc on a notre theta final
            /*
            for(double theta = theta_start ; theta <= theta_end ; theta+=dtheta)    // Qaund on mettra le trou noir, il faudra mettre ici le calcul de r qui sera en fonction de theta
            {}
            */

            double r = b/sin(theta_end);                                    // On calcule r qui dépent du dernier theta, lorsqu'on touche le disque, ici, b/sin fonctionne parce qu'on n'a pas encore de trou noir
            
            if(r>R_min && r<R_max){                                         // Condition pour déterminer si on se trouve bel et bien dans la région du disque d'accrétion
                fprintf(file,"%.20e \t\t %.20e \t\t %.20e\n", b, alpha, r); // si c'est le cas, on écrit la valeur dans le fichier
            }
        }

    }
    printf("\nTerminé !\n");    // Le calcul est terminé
    fclose(file);               // On ferme le fichier à la fin de notre code
}



////////////////////////////////////////////////////////////////////// MAIN //////////////////////////////////////////////////////////////////////
int main(){
    FILE* file = fopen("image_detecteur.txt","w");  // On ouvre un fichier .txt pour écrire dedans

    const double rs = 1;                            // Rayon de Schwarzschild

    const double i = M_PI/6;                        // Angle i, inclinaison du disque d'accrétion
    const double b_start = 0.01;                    // Première valeur de b du détecteur
    const double b_end = 10*rs ;                    // Dernière valeur de b (on arrête b à cette valeur), c'est la taille du détecteur
    const double db = 1e-2;                         // Le pas de b, ce qui va déterminer le nombre d'élément b, + il est petit, + les calculs sont longs

    const double alpha_start = 0;                   // Angle alpha du détecteur, première valeur de alpha
    const double dalpha = 1e-2;                     // Le pas de alpha, idem que db
    
    //const double theta_start = 1e-2;              // On ne l'utilise pas ici, mais il servira pour le disque+trou noir, c'est la première valeur de l'angle theta 
    //const double dtheta = 1e-2;                   // On ne l'utilise pas ici, mais il servira pour le disque+trou noir, c'est le pas de theta, idem que db et dalpha
    
    const double R_min = rs;                        // Rayon du disque d'accrétion minimum qu'on fixe à rs pour avoir trou noir -> disque autour du trou noir
    const double R_max = 7*rs;                      // Rayon maximal du disque d'accrétion

    explorationCI(file, i, b_start, b_end, db, alpha_start, dalpha, R_min, R_max);  // On appelle notre fonction pour faire les calculs

    return 0;                                       // On return 0 le main comme d'hab

}