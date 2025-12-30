#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <omp.h>    // Permet le calcul parallèle



//////////////////////////////////////////////// Fonction qui permet d'afficher une barre de progression ////////////////////////////////////////////////
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



//////////////////////////////////////////////////////////////// Equation des géodésiques ////////////////////////////////////////////////////////////////
double du_point_dtheta(double u) {  // u'' = 3/2 u^2 - u 
	return ((3.0/2.0) * u * u) - u; // Notre équation de trajectoire
}



////////////////////////////////////////////////////////////////// Runge Kutta ordre 4 //////////////////////////////////////////////////////////////////
double RK4(double *u, double *u_point, double dtheta, double rs, double *y, double distance_bg, double *theta) {
    double r=rs/ *u;                                        // On calcule ce qui nous intéresse pour le futur plot

	while (r*cos(*theta)>=distance_bg) {                            // Boucle jusqu'à dépasser le mur invisible derrière le trou noir
		double k1 = (*u_point) * dtheta;                            // Premier terme u de RK4
		double k1_point = du_point_dtheta(*u) * dtheta;             // Premier terme u' de RK4

		double k2 = (*u_point + k1_point/2.0) * dtheta;             // Deuxième terme u de RK4
		double k2_point = du_point_dtheta(*u + k1/2.0) * dtheta;    // Deuxième terme u' de RK4

		double k3 = (*u_point + k2_point/2.0) * dtheta;             // Troisième terme u de RK4
		double k3_point = du_point_dtheta(*u + k2/2.0) * dtheta;    // Troisième terme u' de RK4

		double k4 = (*u_point + k3_point) * dtheta;                 // Quatrième terme u de RK4
		double k4_point = du_point_dtheta(*u + k3) * dtheta;        // Quatrième terme u' de RK4
		
		*u = *u + (k1 + 2*k2 + 2*k3 + k4)/6.0;                                        // u final de RK4
		*u_point = *u_point + (k1_point + 2*k2_point + 2*k3_point + k4_point)/6.0;    // u' final de RK4

        r=rs/ *u;                               // On calcule r à la fin de notre tour de boucle     
        *theta += dtheta;                       // On calcule theta à chaque fois
        *y = r*sin(*theta);                     // On calcule le Y local sur le mur du fond, obligatoire pour retracer les bons pixels avec python
        if (r<rs || *theta>3.0*M_PI/2.0) break; // On regarde si on est tombé dans le trou noir ou si le rayon à fait un tour du trou noir (quasi pas visible sur l'image finale)
        
    }
    return r;   // On retourne r pour l'exploiter en dehors de la fonctoin
}



////////////////////////////////////////////////////////////////////// Fonction Exploration des Conditions Initiales //////////////////////////////////////////////////////////////////////
void explorationCI(double b_end, double b_start, double db, double alpha_start, double daplha, double theta_start, double dtheta, double rs, double distance_bg) {

    FILE* file = fopen("image_detecteur.bin","wb");     // On ouvre un fichier binaire pour que ce soit léger
    if (!file) { perror("fopen"); exit(1); }            // On regarde si le fichier peu s'ouvrir


    const double bc = 3.0 * sqrt(3.0) * rs /2.0;        // constante sans dimension donnée dans le sujet
    
    int n_steps = (int)((b_end - b_start) / db);        // Pour la barre de progression, calcule le nombre total 
    
    int progress_counter = 0;                           // Compteur pour la barre
    
    #pragma omp parallel                                // Début de la région parallèle : les threads OpenMP sont lancés ici
    {
        #pragma omp for schedule(dynamic)                       // La boucle sur 'step' est répartie entre les threads, avec un ordonnancement dynamique
        for (int step = 0; step < n_steps; ++step) {                    // step est privé à chaque thread, géré par OpenMP
            double b = b_start + step * db;                             // On reconstruit b à partir de step
            
            #pragma omp atomic                                          // Incrément atomique: évite les conditions de course sur 'progress_counter'
            progress_counter++;
            
            #pragma omp critical(progress_print)                        // Affichage de la barre pour une seule thread à la fois
            {
                progress_bar(progress_counter / (double)n_steps);                 // Barre de progression qui se met à jour
            }
            
            for(double alpha = alpha_start ; alpha<M_PI ; alpha+=daplha){                               // Boucle sur alpha (quart d'image), (séquentielle à l'intérieur de chaque thread)                     
                double uc = bc/b;                                                                               // Variable sans dimension donnée dans le sujet
                double u = 1e-3;  						    				                                    // r0 = rs/u0; lumière qui arrive de loin, r0 -> inf, donc u0 << 1 ici on prend par ex 1e-3
                double u_point = sqrt( (4.0/27.0) * uc*uc - u*u + u*u*u);                                       // Première valeur de u', à partir de l'équation adimensionnée
                double y=0;                                                                                     // On initialise y à 0
                double theta = theta_start;                                                                     // On initialise theta à notre valeur de départ (normalement 0.0)

                double r = RK4(&u, &u_point, dtheta, rs, &y, distance_bg, &theta);                              // On récupère la dernière valeur de r après avoir résolu l'equa diff 
                
                if(r>rs && theta<3.0*M_PI/2.0){                                                         		    // Condition pour savoir si on n'est pas dans le trou noir, et si qu'on n'a pas fait un tour de trou noir (peu visible)

                    #pragma omp critical(file_write)                                                                // Section critique: un seul thread à la fois peut écrire dans le fichier
                    {
                        float rec[3];                                                                               // On initialise un tableau de taille 3
                        
                        // 1er quart d'image
                        rec[0] = (float)b;                                                                          // On met b en float pour gagner de l'espace, et on le met dans la premiere valeur du tableau
                        rec[1] = (float)alpha;                                                                      // Pareil avec alpha
                        rec[2] = (float)y;                                                                          // Pareil avec y
                        fwrite(rec, sizeof(float), 3, file);                                                        // On écrit dans le fichier notre tableau (c'est du binaire donc pas de fprintf)

                        // 2eme quart d'image
                        rec[0] = (float)b;                                                                          
                        rec[1] = (float)(-alpha + M_PI);
                        rec[2] = (float)y;                              // En fait, on a une symétrie sur x mais aussi sur y, car pas de disque d'accrétion ici, seulement une lentille
                        fwrite(rec, sizeof(float), 3, file);            // Pour ça qu'on peut réduire le temps de calcul par 4 cette fois-ci

                        // 3eme quart d'image
                        rec[0] = (float)b;
                        rec[1] = (float)(alpha + M_PI);
                        rec[2] = (float)y;
                        fwrite(rec, sizeof(float), 3, file);

                        // 4eme quart d'image
                        rec[0] = (float)b;
                        rec[1] = (float)(-alpha);
                        rec[2] = (float)y;
                        fwrite(rec, sizeof(float), 3, file);

                    }
                }
            }
        }
    }
    printf("\nTerminé !\n");    // Tous les threads ont fini leurs tâches
    fclose(file);               // On n'oublie pas de fermer le fichier une fois terminé
}



////////////////////////////////////////////////////////////////////// MAIN //////////////////////////////////////////////////////////////////////
int main(){
    double start = omp_get_wtime();         // Début chronométrage (temps réel)

    const double rs = 1;                    // Rayon de Schwarzschild

    const double b_start = 0.01;            // Première valeur de b du détecteur
    const double b_end = 70*rs;             // Dernière valeur de b (on arrête b à cette valeur), c'est la taille du détecteur
    const double db = 1e-2;                 // Le pas de b, ce qui va déterminer le nombre d'élément b, + il est petit, + les calculs sont longs
    
    const double alpha_start = M_PI/2.0;    // Angle alpha du détecteur, première valeur de alpha (on ne fait qu'une demie image)
    const double daplha = 1e-3;             // Le pas de alpha, idem que db
    
    const double theta_start = 0.0;         // Première valeur de l'angle theta
    const double dtheta = 1e-3;             // Le pas de theta, idem que db et dalpha
    
    const double distance_bg = -30;         // la distance du mur derrière le trou noir, va définir un peu ce qu'on voit avec la lentille, dans la réalité on n'a pas un mur mais plein d'étoiles à des distances différentes

    explorationCI(b_end, b_start, db, alpha_start, daplha, theta_start, dtheta, rs, distance_bg);  // Appel de la fonction qui fait les calculs
    
    
    double end = omp_get_wtime();                               // Fin chronométrage (temps réel)
    double elapsed = end - start;                               // Calcul du temps mis (en secondes)
    printf("Temps total : %.3f secondes\n", elapsed);           // Printf du temps total
    #pragma omp parallel
    {
        #pragma omp single                                      // Une seule thread affiche l’information sur le nombre de threads utilisés
        {
            printf("Threads utilisés : %d\n", omp_get_num_threads());
        }
    }

    return 0;                               // On retourne 0 à la fin du main
}
