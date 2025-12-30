import matplotlib.pyplot as plt			# On importe la librairie pour plot
import numpy as np						# On importe numpy pour les maths


#Constantes (SI)
G = 6.67430e-11      					# Constante grav en m^3.kg^-1.s^-2
c = 299792458        					# Vitesse de la lumière en m.s^-1
m_sun = 1.98847e30   					# Masse du soleil en kg

#Paramètres
m = 5*m_sun								# Masse du trou noir
r_s = (2*G*m)/c**2   					# Rayon de Schwarzschild (m)


u_vals = np.linspace(0,2,100) 			# u = (r_s/r) adimensionnement


def f(U,u_c):
    return (U**2)*(U-1) + (u_c**2)*4/27 # Ma fonction de calcul avec des u et u_c différents



liste_couleurs = ["red", "blue", "green"]					# Couleurs de mes plots
liste_labels = [r"$u_c$ < 1", r"$u_c$ = 1", r"$u_c$ > 1"]	# Labels de mes plots
liste_f = [f(u_vals,0.5), f(u_vals,1), f(u_vals,5)]			# Liste de mes fonction avec des valeurs différentes de u_c


plt.figure(figsize=(10, 8))														# Taille de mes plots
for i in range(3):
	
	plt.plot(u_vals, liste_f[i], color=liste_couleurs[i], label=liste_labels[i])	# Je fais un plot f(U, uc) en fonction de U
	plt.axhline(y=0, color='black', linestyle='--', linewidth=0.7)					# Je plot une ligne pour visualiser le rayon de Schwarzschild
	plt.xlabel("u", fontsize=14)													# On nomme l'axe des x
	plt.ylabel(r"f(u, $u_c$)", fontsize=14)											# On nomme l'axe des y
	plt.grid()
	plt.legend(fontsize=14)															# On ajoute la légende des figure
	plt.tight_layout()																# Ajuster l'espacement entre les labels etc...

plt.show()																			# Afficher les 3 figures d'un coup
