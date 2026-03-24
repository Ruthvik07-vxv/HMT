import sys
import numpy as np
import matplotlib.pyplot as plt


##  Creating a mesh of nx, ny  ##
def CreateMesh(nx, ny, tTop = 0, tBottom = 0, tLeft = 0, tRight = 0) :

    tMesh = np.zeros((ny, nx), dtype=float)
    tMesh[0, :] = tBottom
    tMesh[ny-1, :] = tTop
    tMesh[:, 0] = tLeft
    tMesh[:, nx-1] = tRight
    print("Mesh Created Successfully!")

    return tMesh


# Iterator for converging the solution to tolerance or convergance limit ##
def Iterator(tMesh, nx, ny, tolerance, h=0, k=1, Tinf=0, dirichlet_top=True, dirichlet_bottom=True, dirichlet_left=True, dirichlet_right=True) :

    iterations = 0
    maxIterations = 100000
    error = 100

    dx = 1
    Bi = h*dx/k

    while error > tolerance and iterations < maxIterations:

        error = 0
        iterations += 1
        tOldMesh = tMesh.copy()

        # Interior nodes (pure conduction)
        for j in range(1, nx-1):
            for i in range(1, ny-1):

                tMesh[i, j] = 0.25 * (
                    tMesh[i+1, j] +
                    tMesh[i-1, j] +
                    tMesh[i, j+1] +
                    tMesh[i, j-1]
                )

                change = abs(tMesh[i, j] - tOldMesh[i, j])
                error = max(error, change)

        # Bottom boundary Convection
        if not dirichlet_bottom:
            for j in range(1, nx-1):
                tMesh[0, j] = (
                    tMesh[0, j+1] +
                    tMesh[0, j-1] +
                    tMesh[1, j] +
                    Bi*Tinf
                ) / (3 + Bi)

        # Top boundary Convection
        if not dirichlet_top:
            for j in range(1, nx-1):
                tMesh[ny-1, j] = (
                    tMesh[ny-1, j+1] +
                    tMesh[ny-1, j-1] +
                    tMesh[ny-2, j] +
                    Bi*Tinf
                ) / (3 + Bi)

        # Left boundary Convection
        if not dirichlet_left:
            for i in range(1, ny-1):
                tMesh[i, 0] = (
                    tMesh[i+1, 0] +
                    tMesh[i-1, 0] +
                    tMesh[i, 1] +
                    Bi*Tinf
                ) / (3 + Bi)

        # Right boundary Convection
        if not dirichlet_right:
            for i in range(1, ny-1):
                tMesh[i, nx-1] = (
                    tMesh[i+1, nx-1] +
                    tMesh[i-1, nx-1] +
                    tMesh[i, nx-2] +
                    Bi*Tinf
                ) / (3 + Bi)

        # Corners (only if convection present)

        if not dirichlet_bottom and not dirichlet_right:
            tMesh[0, nx-1] = (
                tMesh[0, nx-2] +
                tMesh[1, nx-1] +
                2*Bi*Tinf
            ) / (2 + 2*Bi)

        if not dirichlet_bottom and not dirichlet_left:
            tMesh[0, 0] = (
                tMesh[0, 1] +
                tMesh[1, 0] +
                2*Bi*Tinf
            ) / (2 + 2*Bi)

        if not dirichlet_top and not dirichlet_right:
            tMesh[ny-1, nx-1] = (
                tMesh[ny-1, nx-2] +
                tMesh[ny-2, nx-1] +
                2*Bi*Tinf
            ) / (2 + 2*Bi)

        if not dirichlet_top and not dirichlet_left:
            tMesh[ny-1, 0] = (
                tMesh[ny-1, 1] +
                tMesh[ny-2, 0] +
                2*Bi*Tinf
            ) / (2 + 2*Bi)

        if iterations % 500 == 0:
            print()
            print("Iterations:", iterations)
            print("Error:", error)

    print("Total Iterations taken:", iterations)
    print("Final error:", error)

    return tMesh


def PrintTemperatureGrid(tMesh):
    for row in tMesh:
        for val in row:
            print(f"{val:8.2f}", end=" ")
        print()


## Theoretical Solution (only forr conduction) ##
def theoreticalSolution(x, y, Lx, Ly, tTop, tBottom, tLeft, tRight, terms = 100) :
    # Let t0 = tBottom for the calculation #
    T = 0 
    tFinal_Bottom = tBottom 
    tFinal_Top = 0
    tFinal_Left = 0
    tFinal_Right = 0

    tShift_top = tTop 

    for n in range(1, terms, 2) :

        # Top Layer #
        tFinal_Top += ((4 * (tTop - tBottom) / (n * np.pi)) * 
              np.sinh(n * np.pi * y / Lx) /
              np.sinh(n * np.pi * Ly / Lx) *
              np.sin(n * np.pi * x / Lx))
        
        # Left Layer #
        tFinal_Left += ((4 * (tLeft - tBottom) / (n * np.pi)) *
              np.sinh(n * np.pi * (Lx - x) / Ly) /
              np.sinh(n * np.pi * Lx / Ly) *
              np.sin(n * np.pi * y / Ly))
        
        # Right Layer #
        tFinal_Right += ((4 * (tRight - tBottom) / (n * np.pi)) *
              np.sinh(n * np.pi * x / Ly) /
              np.sinh(n * np.pi * Lx / Ly) *
              np.sin(n * np.pi * y / Ly))
        
    T = tFinal_Top + tFinal_Bottom + tFinal_Left + tFinal_Right
    return T


## Theoretical Solution Grid ##
def analyticalGrid(nx, ny, Lx, Ly, tTop, tBottom, tLeft, tRight) :
    tTheory = np.zeros((ny, nx))
    for i in range(0, ny) :
        for j in range (0, nx) :
            x = j / (nx - 1) * Lx
            y = i / (ny - 1) * Ly

            tTheory[i, j] = theoreticalSolution(x, y, Lx, Ly, tTop, tBottom, tLeft, tRight)
    
    return tTheory


## Plotting temperature contours ##
def PlotTemperature(tMesh, name, plottype) :
    plt.contourf(tMesh, 200, cmap = "inferno")
    plt.colorbar()
    plt.title(plottype)
    plt.savefig(name, dpi=300)
    print("Temperature Countours is saved as \"Temperature Countour.png\"")
    plt.show()


def PlotIsotherms(tMesh, name, plottype) :
    cs = plt.contour(tMesh, 20)
    plt.clabel(cs, inline = True, fontsize = 8)
    plt.title(plottype)
    plt.savefig(name, dpi=300)
    print("Temperature Isotherms is saved as \"Temperature Isotherms.png\"")
    plt.show()


def PlotCombined(tMesh) :
    plt.contourf(tMesh, 80, cmap = "inferno")
    plt.colorbar()
    cs = plt.contour(tMesh, 20, colors="black")
    plt.clabel(cs, inline = True, fontsize = 8)
    plt.title("Temperature Distribution with isotherms")
    plt.savefig("Combined Plot.png", dpi=300)
    print("Combined plot is saved as \"Combined Plot.png\"")
    plt.show()


## Saving temperature grid as .txt file ##
def saveTemperatureGrid(tMesh, filename = "TemperatureGrid.txt") :
    with open(filename, "w") as f :
        for row in np.flipud(tMesh) :
            for val in row :
                f.write(f"{val:10.4f} ")
            f.write("\n")
    
    print(f"Temperature Grid saved to {filename}")


def main(nx = 40, ny = 40, tTop = 100, tBottom = 0, tLeft = 0, tRight = 0, tolerance = 1e-6, Lx = 1, Ly = 1) :
    conv = input("Is there is convection at the boundaries(y/n): ").lower()
    if conv == 'y' :
        h = float(input("Enter the convective coefficient h: "))
        k = float(input("Enter the thermal conductivity of the material k: "))
        tInf = float(input("Enter the ambient temperature: "))
        dirichlet_top = input("Is there is convection at the top surface(y/n): ").lower() == 'n'
        dirichlet_bottom = input("Is there is convection at the Bottom surface(y/n): ").lower() == 'n'
        dirichlet_left = input("Is there is convection at the Left surface(y/n): ").lower() == 'n'
        dirichlet_right = input("Is there is convection at the Right surface(y/n): ").lower() == 'n'
    else :
        h = 0
        k = 1
        tInf = 0
        dirichlet_top = True
        dirichlet_bottom = True
        dirichlet_left = True
        dirichlet_right = True

    tMesh = CreateMesh(nx, ny, tTop, tBottom, tLeft, tRight)
    tMesh = Iterator(tMesh, nx, ny, tolerance, h, k, tInf, dirichlet_top, dirichlet_bottom, dirichlet_left, dirichlet_right)

    PlotTemperature(np.fliplr(tMesh), "Temperature Countours.png", "Temperature Contours")
    PlotIsotherms(np.fliplr(tMesh), "Temperature Isotherms.png", "Temperature Isotherms")
    PlotCombined(np.fliplr(tMesh))
    saveTemperatureGrid(tMesh, 'ActualTemperature.txt')
    np.savetxt("TemperatureField.csv", np.flipud(tMesh), delimiter=",")
    print()
    print()

    if conv != 'y' :
        print("Analytical Solution: ")
        tTheory = analyticalGrid(nx, ny, Lx, Ly, tTop, tBottom, tLeft, tRight)
        print("Theoretical temperatures are created successfully")
        saveTemperatureGrid(tTheory, 'TheoreticalT.txt')
        PlotTemperature(tTheory, "Theoretical temperature.png", "Temperature Contours")
        PlotIsotherms(tTheory, "Theoretical isotherms.png", "Temperature Isotherms")
        np.savetxt("TheoreticalT.csv", np.flipud(tTheory), delimiter=",")
        print("Theoretical Temperatures are saved as TheoreticalT.txt and TheoreticalT.csv respectively")
        errorActual = np.abs(tMesh - tTheory)
        print("Errors")
        print(f"max error: {np.max(errorActual)}")
        print(f"Average error: {np.average(errorActual)}")
        np.savetxt("errorMap.csv", np.flipud(errorActual), delimiter=",")
        saveTemperatureGrid(errorActual, "errorMap.txt")
        PlotTemperature(errorActual, "Error Graph.png", "Error Mapping")
        print("Error plot values is saved as errorMap.csv and errorMap.txt respectively")
        
        print(f"Enter the values of x (< {nx}) and y (< {ny}) to compare the values")
        verify = 1
        while verify != -1 :
            np.allclose(tTheory, tMesh.T)
            x = int(input(f"Enter the value of x (< {nx}): "))
            y = int(input(f"Enter the value of y (< {ny}): "))
            print(f"Analytical Solution: {tTheory[y, x]}")
            print(f"Converged Solution: {tMesh[y, x]}")
            print(f"Error in Solution: {errorActual[y, x]}")
            verify = int(input("Enter any value to continue, and -1 to terminate the process: "))
            if verify == -1 :
                sys.exit()
            else :
                continue
    verify = 1
    while verify != -1 :
        x = int(input(f"Enter the value of x (< {nx}): "))
        y = int(input(f"Enter the value of y (< {ny}): "))
        print(f"Converged Solution: {tMesh[y, x]}")
        verify = int(input("Enter any value to continue, and -1 to terminate the process: "))
        if verify == -1 :
                sys.exit()
        else : 
            continue

Lx = float(input("Enter the length of the specimen: "))
Ly = float(input("Enter the Height of the specimen: "))
print("Define the total number of rows and columns (-1 to auto assign)")
nx = int(input("Enter the total number of rows: "))
if nx == -1 :
    nx = 40
ny = int(input("Enter the total number of columns: "))
if ny == -1 :
    ny = 40
print()
print("Enter the temperatures of the 4 egdes (-1 if not specified)")
tTop = float(input("Enter the temperature of the Top surface: "))
if tTop == -1 :
    tTop = 0
if tTop <= -273.15 :
    print("Error in temperature, cannot be less than -273.15 in celsius")
    sys.exit()
tBottom = float(input("Enter the temperature of the Bottom surface: "))
if tBottom == -1 :
    tBottom = 0
if tBottom <= -273.15 :
    print("Error in temperature, cannot be less than -273.15 in celsius")
    sys.exit()
tLeft = float(input("Enter the temperature of the Left surface: "))
if tLeft == -1 :
    tLeft = 0
if tLeft <= -273.15 :
    print("Error in temperature, cannot be less than -273.15 in celsius")
    sys.exit()
tRight = float(input("Enter the temperature of the Right surface: "))
if tRight == -1 :
    tRight = 0
if tRight <= -273.15 :
    print("Error in temperature, cannot be less than -273.15 in celsius")
    sys.exit()
tolerance = float(input("Enter the Convergence Error: "))


main(nx, ny, tTop, tBottom, tLeft, tRight, tolerance, Lx, Ly)
