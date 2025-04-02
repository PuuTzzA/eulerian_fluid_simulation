const canvas = document.getElementById("canvas");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
const ctx = canvas.getContext("2d");

window.addEventListener("resize", e => {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
})

class GraphicsSettings {
    constructor() {
        this.lineWidth = 2;
        this.lineColour = "rgb(255, 255, 255)";
        this.particleSize = 3;
        this.minSpeed = 0;
        this.maxSpeed = 2000;
        this.velocityRange = [350, 110]; // "colorramp" [hue_for_min_speed, hue_for_max_speed]
    }

    getLineWidth() {
        return this.lineWidth;
    }

    getLineColour() {
        return this.lineColour;
    }

    getParticleRadius() {
        return this.particleSize;
    }

    getParticleColour(cell) {
        let hue = 300;
        let particle = {selected: false};
        return particle.selected ? "white" : `hsla(${hue}, 100%, 50%, 1)`;
    }
}

class Fluid {
    constructor(resolutionX, resolutionY) {
        this.gs = new GraphicsSettings();

        this.resolutionX = resolutionX;
        this.resolutionY = resolutionY;

        this.mousePos = [0, 0];
        this.mousePosPrev = [0, 0];
        this.mouseRadius = 100;
        this.mouseRadiusSquared = this.mouseRadius * this.mouseRadius;
        this.mousePressed = false;
        this.mouseForce = 200;

        canvas.width = window.innerWidth;
        canvas.height = window.innerHeight;

        this.grid;
        this.gridSize = 10;
    }

    setResolutionX(newX) {
        this.initGrid(newX, this.resolutionY);
    }

    setResolutionY(newY) {
        this.initGrid(this.resolutionX, newY);
    }

    initGrid(resX, resY) {
        let shouldX =  window.innerWidth / resX;
        let shouldY =  window.innerHeight / resY;

        this.gridSize = Math.min(shouldX, shouldY);
        this.resolutionX = resX;
        this.resolutionY = resY;
    }

    update(dt) {
    }

    draw() {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        ctx.strokeStyle = this.gs.getLineColour();
        ctx.lineWidth = this.gs.getLineWidth();

        const offsetTop = window.innerHeight - this.gridSize * this.resolutionY;

        for (let x = 0; x < this.resolutionX; x++) {
            for (let y = 0; y < this.resolutionY; y++) {

                ctx.fillStyle = this.gs.getParticleColour(10);
    
                ctx.strokeRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);
                ctx.fillRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);
            }
        }
    }
}
