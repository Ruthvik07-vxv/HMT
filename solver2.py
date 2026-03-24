import numpy as np
import matplotlib.pyplot as plt

def CreateMesh(nx, ny, tTop = 0, tBottom = 0, tLeft = 0, tRight = 0) :
    ##  Creating a mesh of nx, ny  ##
    tMesh = np.zeros((ny, nx), dtype=float)

    ## Dirichlet Conditions ##
    tMesh[0, :] = tBottom
    tMesh[ny-1, :] = tTop
    tMesh[:, 0] = tLeft
    tMesh[:, nx-1] = tRight

    ##  Creating a fixed temperature Node  ##
    fixed = np.zeros((ny, nx), dtype=bool)

    fixed[0, :] = True
    fixed[ny-1, :] = True
    fixed[:, 0] = True
    fixed[:, nx-1] = True
    print("Mesh Created Successfully!")

    return tMesh, fixed

def fixedTemperatures(tMesh, fixed, x, y, nx, ny, temp) :
    if x > nx-1 or x < 0 :
        print("Error! value of X coordinate is greater than the number of cells")
        x = int(input("Enter the value of x again: "))
        return fixedTemperatures(tMesh, fixed, x, y, nx, ny, temp)
    elif y > ny-1 or y < 0 :
        print("Error! value of Y coordinate is greater than the number of cells")
        y = int(input("Enter the value of y again: "))
        return fixedTemperatures(tMesh, fixed, x, y, nx, ny, temp)
    else :
        tMesh[y, x] = temp
        fixed[y, x] = True
        print("Value added Successfully")

    return tMesh, fixed

def Iterator_NoHeatGeneration (tMesh, fixed, nx, ny, tolerance, error = 100) :
    iterations = 0
    maxIterations = 100000
    while error > tolerance and iterations < maxIterations :
        ##  To find the new max error  ##
        error = 0 
        iterations += 1
        ##  Copy the old data  ##
        tOldMesh = tMesh.copy()

        ##  Iterate until convergence  ##
        for j in range(1, nx-1) :
            for i in range(1, ny-1) :
                if fixed[i, j] :
                    continue
                    
                tMesh[i, j] = 0.25 * (tMesh[i+1, j] + tMesh[i-1, j] + tMesh[i, j+1] + tMesh[i, j-1])
                change = abs(tMesh[i, j] - tOldMesh[i, j])
                error = max(error, change)

        if iterations % 500 == 0 :
            print()
            print("Iterations: ", iterations)
            print("Error: ", error)

        if iterations == maxIterations :
            print("Max Iterations reached! Please try a lower Convergence value")
            print("Final error without convergence: ", error)
            print("Convergence Value: ", tolerance)
            return tMesh

    
    print("Total Iterations taken: ", iterations)
    print("Final error after Convergence: ", error)
    return tMesh

def Iterator_HeatGeneration (tMesh, fixed, nx, ny, tolerance, length, height, q, k, error = 100) :
    dx = length / (nx - 1)
    dy = height / (ny - 1)
    iterations = 0
    maxIterations = 100000
    while error > tolerance and iterations < maxIterations :
        ##  To find the new max error  ##
        error = 0 
        iterations += 1
        ##  Copy the old data  ##
        tOldMesh = tMesh.copy()

        ##  Iterate until convergence  ##
        for j in range(1, nx-1) :
            for i in range(1, ny-1) :
                if fixed[i, j] :
                    continue
                    
                tMesh[i, j] = (((tMesh[i+1, j] + tMesh[i-1, j]) / dx**2) + 
                               ((tMesh[i, j+1] + tMesh[i, j-1]) / dy**2) + 
                               (q/k)) / ((2 / dx**2) + (2 / dy**2))
                change = abs(tMesh[i, j] - tOldMesh[i, j])
                error = max(error, change)
        
        if iterations % 500 == 0 :
            print("Iterations: ", iterations)
            print("Error: ", error)
            print()

        if iterations == maxIterations :
            print("Max Iterations reached! Please try a lower Convergence value")
            print("Final error without convergence: ", error)
            print("Convergence Value: ", tolerance)
            return tMesh

    print("Total Iterations taken: ", iterations)
    print("Final error after Convergence: ", error)
    return tMesh

def PlotTemperature(tMesh) :
    plt.contourf(tMesh, 200, cmap = "inferno")
    plt.colorbar()
    plt.title("Temperature Field")
    plt.savefig("Temperature Contours.png", dpi=300)
    print("Temperature Countours is saved as \"Temperature Countour.png\"")
    plt.show()


def PlotIsotherms(tMesh) :
    cs = plt.contour(tMesh, 20)
    plt.clabel(cs, inline = True, fontsize = 8)
    plt.title("Temperature Isotherms")
    plt.savefig("Temperature Isotherms.png", dpi=300)
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



def main_NoHeatGeneration(nx, ny, tTop, tBottom, tLeft, tRight, tolerance) :
    tMesh, fixed = CreateMesh(nx, ny, tTop, tBottom, tLeft, tRight)
    print("Enter Y or y if there are any fixed temperatures in the mesh ", nx - 1, " ", ny - 1, ". If there are none present, click any other key")
    test = input("Enter a key: ")
    if test == "Y" or test == "y" :
        while (1) :
            x = int(input("Enter the value of x: "))
            y = int(input("Enter the value of y: "))
            temp = int(input("Enter the temperature to the node: "))

            fixedTemperatures(tMesh, fixed, x, y, nx, ny, temp)
            print("Enter any value to continue the operation, and -1 to complete the operation i.e, when all fixed temperatures are given")
            cnt = int(input("Enter a value (-1 if operation ended): "))
            if cnt == -1 :
                break
            else :
                continue
    tMesh = Iterator_NoHeatGeneration(tMesh, fixed, nx, ny, tolerance) 
    PlotTemperature(tMesh)
    PlotIsotherms(tMesh)
    PlotCombined(tMesh)

def main_HeatGeneration (nx, ny, tTop, tBottom, tLeft, tRight, tolerance, q, k, Length, Height) :
    tMesh, fixed = CreateMesh(nx, ny, tTop, tBottom, tLeft, tRight)
    print("Enter Y or y if there are any fixed temperatures in the mesh ", nx - 1, " ", ny - 1, ". If there are none present, click any other key")
    test = input("Enter a key: ")
    if test == "Y" or test == "y" :
        while (1) :
            x = int(input("Enter the value of x: "))
            y = int(input("Enter the value of y: "))
            temp = int(input("Enter the temperature to the node: "))

            fixedTemperatures(tMesh, fixed, x, y, nx, ny, temp)
            print("Enter any value to continue the operation, and -1 to complete the operation i.e, when all fixed temperatures are given")
            cnt = int(input("Enter a value (-1 if operation ended): "))
            if cnt == -1 :
                break
            else :
                continue
    tMesh = Iterator_HeatGeneration(tMesh, fixed, nx, ny, tolerance, Length, Height, q, k) 
    PlotTemperature(tMesh)
    PlotIsotherms(tMesh)
    PlotCombined(tMesh)

print("Enter the type of the heat transfer thats taking place")
print("Type 1 if the heat transfer that's taking place is Steady state heat conduction with no heat generation")
print("Type 2 if the heat transfer that's taking place is Steady state heat conduction with heat generation")
HT = int(input("Enter the case: "))

if HT == 1 :
    nx = int(input("Enter the total number of cells used on x axis: ")) + 1
    ny = int(input("Enter the total number of cells used on y axis: ")) + 1
    tTop = float(input("Enter the temperature of the Top surface: "))
    tBottom = float(input("Enter the temperature of the Bottom surface: "))
    tLeft = float(input("Enter the temperature of the Left surface: "))
    tRight = float(input("Enter the temperature of the Right surface: "))
    tolerance = float(input("Enter the Convergence Error: "))
    main_NoHeatGeneration(nx, ny, tTop, tBottom, tLeft, tRight, tolerance)

elif HT == 2 :
    nx = int(input("Enter the total number of cells used on x axis: ")) + 1
    ny = int(input("Enter the total number of cells used on y axis: ")) + 1
    tTop = float(input("Enter the temperature of the Top surface: "))
    tBottom = float(input("Enter the temperature of the Bottom surface: "))
    tLeft = float(input("Enter the temperature of the Left surface: "))
    tRight = float(input("Enter the temperature of the Right surface: "))
    Length = float(input("Enter the length of the object used: "))
    Height = float(input("Enter the height of the object used: "))
    q = float(input("Enter the amount of the heat being generated in the objece: "))
    k = float(input("Enter the Thermal Conductivity Constant of the Material: "))
    tolerance = float(input("Enter the Convergence Error: "))
    main_HeatGeneration(nx, ny, tTop, tBottom, tLeft, tRight, tolerance, q, k, Length, Height)

else :
    print("Wrong Choice! Try again")