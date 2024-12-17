//SIMULATION PARAMETERS: Change these values to change the simulation settings
//Note that all variables in this simulation are stored as 64-bit floating-point values, meaning that all variables have 16 significant digits in the backgroud of the simulation, regardless of the significance of the input variables

//Thermal conductivity of medium (in W m^-1 K-1)
thermConductivity = 6.2;

// Length step (voxel edge length) (1 = 1m, 0.001=1mm)
steplength = 0.01;

//time step (amount of modeled time that passes between each "frame" of the simulation) in seconds
t = 0.5;

// Size of reactor (width, height) multiples of the legnth step, or units (refered to as voxels)
// Notes that depth is fixed to a singular voxel to best model a 2-dimensional reactor
reactorWidth = 100;
reactorHeight = 100;

// volume of voxel is the length step cubed (m^3) times 1000 (dm^3)
volume = steplength * steplength * steplength * 1000;

//specific heat capacity of medium (in J g-1 K-1)
specheatcapacity = 0.5;

//density of medium per 1dm^3 in g dm^-3
density = 4000;

//specific heat capacity of voxel (the number of Joules of energy it takes to raise the temperature of a voxel by 1 Kelvin)
voxelheat = volume * density * specheatcapacity;

//diffusion constant of reactants and products (in m^2 s^-1)
diffusion = 1 * 10 ** -4;

//enthaply of reaction (in J mol^-1)
dH = -849000;

// activation energy of reaction (in J mol^-1)
Ea = 80000;

//Arrhenius factor

a = 0.001;

//inital temperature of reactor (in K)
initTemp = 298;

//initial tmeperature of "spark"
sparktemp = 9000;

//inital concentration of reactants (in mol dm^3)
initconc = 100;

//display simulation?
displaysim = true;

//Reaction "completion" threshold (in amount of product made relative to inital amount of reactant)
thresh = 0.9;
//Internal variables, do not modify
let rttwo;
thermfluxpre =
  (steplength * steplength * (thermConductivity / steplength) * t) / voxelheat;
thermfluxprert =
  (steplength * steplength * (thermConductivity / (rttwo * steplength)) * t) /
  voxelheat;

difffluxpre = steplength * steplength * (diffusion / steplength) * t;
difffluxprert =
  steplength * steplength * (diffusion / (rttwo * steplength)) * t;

let reactor;
let tempreactor;
let time = 0;
let done = false;
t1 = false;
t2 = false;
t3 = false;
t4 = false;
t5 = false;

class voxel {
  constructor(conc, temp) {
    this.conc = conc;
    this.temp = temp;
    this.pro = 0;
    this.k = 0;
  }

  react() {
    this.k = a * exp((-1 * Ea) / (this.temp * 8.31));
    let concdiff = this.conc * this.conc * this.k * t;
    if (concdiff > this.conc) {
      concdiff = this.conc;
    }
    this.conc -= concdiff;
    this.pro += concdiff;
    let molmade = concdiff * volume;
    this.temp += (-1 * dH * molmade) / voxelheat;
  }
}

function makeReactor(cols, rows, conc, temp) {
  let arr = new Array(cols);
  for (let i = 0; i < arr.length; i++) {
    arr[i] = new Array(rows);
    // Fill the array with class objects
    for (let j = 0; j < arr[i].length; j++) {
      arr[i][j] = new voxel(conc, temp);
    }
  }
  return arr;
}

function setup() {
  rttwo = sqrt(2);
  thermfluxpre = steplength * steplength * (thermConductivity / steplength) * t;
  thermfluxprert =
    steplength * steplength * (thermConductivity / (rttwo * steplength)) * t;

  difffluxpre = steplength * steplength * (diffusion / steplength) * t;
  difffluxprert =
    steplength * steplength * (diffusion / (rttwo * steplength)) * t;

  createCanvas(5 * reactorWidth, 5 * reactorHeight);

  reactor = makeReactor(reactorWidth, reactorHeight, initconc, initTemp);
  tempreactor = Array.from(reactor);

  reactor[0][0].temp = sparktemp;
}

// Check if a coordinate is within the bounds
function withinGrid(x, y) {
  return x >= 0 && x <= reactorWidth - 1 && y >= 0 && y <= reactorHeight - 1;
}

function tempConduct(x, y) {
  let temp = reactor[x][y].temp;

  if (withinGrid(x - 1, y - 1)) {
    temp += thermfluxprert * (reactor[x - 1][y - 1].temp - reactor[x][y].temp);
  }
  if (withinGrid(x + 1, y - 1)) {
    temp += thermfluxprert * (reactor[x + 1][y - 1].temp - reactor[x][y].temp);
  }
  if (withinGrid(x - 1, y + 1)) {
    temp += thermfluxprert * (reactor[x - 1][y + 1].temp - reactor[x][y].temp);
  }
  if (withinGrid(x + 1, y + 1)) {
    temp += thermfluxprert * (reactor[x + 1][y + 1].temp - reactor[x][y].temp);
  }

  if (withinGrid(x - 1, y)) {
    temp += thermfluxpre * (reactor[x - 1][y].temp - reactor[x][y].temp);
  }
  if (withinGrid(x + 1, y)) {
    temp += thermfluxpre * (reactor[x + 1][y].temp - reactor[x][y].temp);
  }
  if (withinGrid(x, y + 1)) {
    temp += thermfluxpre * (reactor[x][y + 1].temp - reactor[x][y].temp);
  }
  if (withinGrid(x, y - 1)) {
    temp += thermfluxpre * (reactor[x][y - 1].temp - reactor[x][y].temp);
  }
  return temp;
}

function reactDiffuse(x, y) {
  temp = reactor[x][y].conc;

  if (withinGrid(x - 1, y - 1)) {
    temp += difffluxprert * (reactor[x - 1][y - 1].conc - reactor[x][y].conc);
  }
  if (withinGrid(x + 1, y - 1)) {
    temp += difffluxprert * (reactor[x + 1][y - 1].conc - reactor[x][y].conc);
  }
  if (withinGrid(x - 1, y + 1)) {
    temp += difffluxprert * (reactor[x - 1][y + 1].conc - reactor[x][y].conc);
  }
  if (withinGrid(x + 1, y + 1)) {
    temp += difffluxprert * (reactor[x + 1][y + 1].conc - reactor[x][y].conc);
  }

  if (withinGrid(x - 1, y)) {
    temp += difffluxpre * (reactor[x - 1][y].conc - reactor[x][y].conc);
  }
  if (withinGrid(x + 1, y)) {
    temp += difffluxpre * (reactor[x + 1][y].conc - reactor[x][y].conc);
  }
  if (withinGrid(x, y + 1)) {
    temp += difffluxpre * (reactor[x][y + 1].conc - reactor[x][y].conc);
  }
  if (withinGrid(x, y - 1)) {
    temp += difffluxpre * (reactor[x][y - 1].conc - reactor[x][y].conc);
  }
  return temp;
}

function proDiffuse(x, y) {
  temp = reactor[x][y].pro;

  if (withinGrid(x - 1, y - 1)) {
    temp += difffluxprert * (reactor[x - 1][y - 1].pro - reactor[x][y].pro);
  }
  if (withinGrid(x + 1, y - 1)) {
    temp += difffluxprert * (reactor[x + 1][y - 1].pro - reactor[x][y].pro);
  }
  if (withinGrid(x - 1, y + 1)) {
    temp += difffluxprert * (reactor[x - 1][y + 1].pro - reactor[x][y].pro);
  }
  if (withinGrid(x + 1, y + 1)) {
    temp += difffluxprert * (reactor[x + 1][y + 1].pro - reactor[x][y].pro);
  }

  if (withinGrid(x - 1, y)) {
    temp += difffluxpre * (reactor[x - 1][y].pro - reactor[x][y].pro);
  }
  if (withinGrid(x + 1, y)) {
    temp += difffluxpre * (reactor[x + 1][y].pro - reactor[x][y].pro);
  }
  if (withinGrid(x, y + 1)) {
    temp += difffluxpre * (reactor[x][y + 1].pro - reactor[x][y].pro);
  }
  if (withinGrid(x, y - 1)) {
    temp += difffluxpre * (reactor[x][y - 1].pro - reactor[x][y].pro);
  }
  return temp;
}

function draw() {
  background(0);
  noStroke();
  for (let x = 0; x <= reactorWidth - 1; x++) {
    for (let y = 0; y <= reactorHeight - 1; y++) {
      reactor[x][y].react();
    }
  }

  for (let x = 0; x <= reactorWidth - 1; x++) {
    for (let y = 0; y <= reactorHeight - 1; y++) {
      tempreactor[x][y].temp = tempConduct(x, y);
    }
  }

  reactor = Array.from(tempreactor);

  if (displaysim) {
    for (let x = 0; x <= reactorWidth - 1; x++) {
      for (let y = 0; y <= reactorHeight - 1; y++) {
        fill(reactor[x][y].temp / 20, reactor[x][y].pro * 5, reactor[x][y].pro);
        square(x * 5, y * 5, 5);
      }
    }
    time += t;
    if (reactor[17][69].pro >= initconc * thresh && !t1) {
      print("P1 speed is: " + (sqrt(17 * 17 + 69 * 69) * steplength) / time);
      t1 = true;
    }
    if (reactor[32][63].pro >= initconc * thresh && !t2) {
      print("P2 speed is: " + (sqrt(32 * 32 + 63 * 63) * steplength) / time);
      t2 = true;
    }
    if (reactor[50][50].pro >= initconc * thresh && !t3) {
      print("P3 speed is: " + (sqrt(50 * 50 * 2) * steplength) / time);
      t3 = true;
    }
    if (reactor[63][32].pro >= initconc * thresh && !t4) {
      print("P4 speed is: " + (sqrt(32 * 32 + 63 * 63) * steplength) / time);
      t4 = true;
    }
    if (reactor[69][17].pro >= initconc * thresh && !t5) {
      print("P5 speed is: " + (sqrt(17 * 17 + 69 * 69) * steplength) / time);
      t5 = true;
    }
    if (reactor[reactorWidth - 1][reactorHeight - 1].pro >= initconc * 0.5) {
      noLoop();
    }
  }
}
