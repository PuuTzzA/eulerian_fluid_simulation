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
        this.maxSpeed = 10;
        this.maxVelocitySize = 1;
        this.arrowWidth = 2;
        this.arrowHeadSize = 9;
        this.arrowColor = "red";
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

    getCellColour(red, green, blue) {
        return `rgb(${red}, ${green}, ${blue})`;
    }

    getArrowWidth() {
        return this.arrowWidth;
    }

    getArrowHeadSize() {
        return this.arrowHeadSize;
    }

    getArrowColor(velocity) {
        return this.arrowColor;
    }

    getArrowLength(velocity) {
        // returns a fraction of the grid spacing
        return (velocity / this.maxSpeed) * this.maxVelocitySize;
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

        this.pressures = new Float32Array(1);
        this.velocitiesX = new Float32Array(1);
        this.velocitiesY = new Float32Array(1);
        this.colorRed = new Float32Array(1);

        this.gridSize = 10;
    }

    setResolutionX(newX) {
        this.initGrid(parseInt(newX), this.resolutionY);
    }

    setResolutionY(newY) {
        this.initGrid(this.resolutionX, parseInt(newY));
    }

    initGrid(resX, resY) {
        let shouldX = window.innerWidth / resX;
        let shouldY = window.innerHeight / resY;

        this.gridSize = Math.min(shouldX, shouldY);
        this.resolutionX = resX;
        this.resolutionY = resY;

        this.pressures = new Float32Array(resX * resY);
        this.colorRed = new Float32Array(resX * resY);
        this.velocitiesX = new Float32Array((resX + 1) * resY);
        this.velocitiesY = new Float32Array(resX * (resY + 1));

        this.velocitiesX[(resX + 1) * 4 + 3] = 2;
        this.velocitiesX[(resX + 1) * 5 + 3] = 2;
        this.velocitiesX[(resX + 1) * 6 + 3] = 2;
    }

    getVelocityAtPoint(x, y) {
        const u = Fluid.bilinearInterpolation(x + 0.5, y, this.velocitiesX, this.resolutionX + 1, this.resolutionY);
        const v = Fluid.bilinearInterpolation(x, y + 0.5, this.velocitiesY, this.resolutionX, this.resolutionY + 1);
        return [u, v];
    }
    
    static bilinearInterpolation(x, y, values, resolutionX, resolutionY) {
        // Clamp to avoid out-of-bounds access
        x = Math.max(0, Math.min(x, resolutionX - 1));
        y = Math.max(0, Math.min(y, resolutionY - 1));

        const x0 = Math.floor(x);
        const y0 = Math.floor(y);
        const x1 = Math.min(x0 + 1, resolutionX - 1);
        const y1 = Math.min(y0 + 1, resolutionY - 1);

        const alpha = x - x0;
        const beta = y - y0;

        const topLeft = values[y0 * resolutionX + x0];
        const topRight = values[y0 * resolutionX + x1];
        const bottomLeft = values[y1 * resolutionX + x0];
        const bottomRight = values[y1 * resolutionX + x1];

        const top = topLeft * (1 - alpha) + topRight * alpha;
        const bottom = bottomLeft * (1 - alpha) + bottomRight * alpha;

        return top * (1 - beta) + bottom * beta;
    }

    static bicubicInterpolation(x, y, values, resolutionX, resolutionY) {

    }

    advect(dt, oldVals, resolutionX, resolutionY) {
        let newVals = new Float32Array(resolutionX * resolutionY);

        for (let i = 0, n = newVals.length; i < n; i++) {
            // if (i != 0) { continue };
            // Runge Kutta 4
            const x = i % resolutionX;
            const y = Math.floor(i / resolutionX);

            const k1 = this.getVelocityAtPoint(x, y);
            const k2 = this.getVelocityAtPoint(x - k1[0] * dt / 2, y - k1[1] * dt / 2);
            const k3 = this.getVelocityAtPoint(x - k2[0] * dt / 2, y - k2[1] * dt / 2);
            const k4 = this.getVelocityAtPoint(x - k3[0] * dt, y - k3[1] * dt);

            let newX = x - (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
            let newY = y - (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

            newX = Math.max(newX, 0);
            newX = Math.min(newX, resolutionX - 1);
            newY = Math.max(newY, 0);
            newY = Math.min(newY, resolutionY - 1);

            newVals[i] = Fluid.bilinearInterpolation(newX, newY, oldVals, resolutionX, resolutionY);
        }

        return newVals;
    }

    addVector(x, y) {
        return [x[0] + y[0], x[1] + y[1]];
    }

    suptractVector(x, y) {
        return [x[0] - y[0], x[1] - y[1]];
    }

    multiplyVector(x, y) {
        return [x[0] * y[0], x[1] * y[1]];
    }

    multiplyVectorScalar(x, a) {
        return [x[0] * a, x[1] * a];
    }

    update(dt) {
        this.colorRed = this.advect(dt, this.colorRed, this.resolutionX, this.resolutionY);
        this.velocitiesX = this.advect(dt, this.velocitiesX, this.resolutionX + 1, this.resolutionY);
        this.velocitiesY = this.advect(dt, this.velocitiesY, this.resolutionX, this.resolutionY + 1);
    }

    draw() {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        ctx.strokeStyle = this.gs.getLineColour();
        ctx.lineWidth = this.gs.getLineWidth();

        const offsetTop = window.innerHeight - this.gridSize * this.resolutionY;


        for (let i = 0, n = this.pressures.length; i < n; i++) {
            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            ctx.fillStyle = this.gs.getCellColour(this.colorRed[i], 0, 0);

            ctx.strokeRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);
            ctx.fillRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);
        }

        let size = 0.4;

        // horizontal velocities
        for (let i = 0, n = this.velocitiesX.length; i < n; i++) {
            if (this.velocitiesX[i] == 0) continue;

            const x = i % (this.resolutionX + 1);
            const y = Math.floor(i / (this.resolutionX + 1));

            let fromX = x * this.gridSize;
            let fromY = (y + 0.5) * this.gridSize + offsetTop;
            let toX = (x + this.gs.getArrowLength(this.velocitiesX[i])) * this.gridSize;
            let toY = fromY;

            let col = i == 0 ? "white" : this.gs.getArrowColor();
            Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
        }

        // vertical velocities
        for (let i = 0, n = this.velocitiesY.length; i < n; i++) {
            if (this.velocitiesY[i] == 0) continue;

            const x = i % (this.resolutionX);
            const y = Math.floor(i / (this.resolutionX));

            let fromX = (x + 0.5) * this.gridSize;
            let fromY = y * this.gridSize + offsetTop;
            let toX = fromX;
            let toY = (y + this.gs.getArrowLength(this.velocitiesY[i])) * this.gridSize + offsetTop;

            let col = i == 0 ? "white" : this.gs.getArrowColor();

            Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
        }
    }

    static drawArrow(ctx, fromx, fromy, tox, toy, arrowWidth, headlen, color) {
        //variables to be used when creating the arrow
        var angle = Math.atan2(toy - fromy, tox - fromx);

        ctx.save();
        ctx.strokeStyle = color;

        //starting path of the arrow from the start square to the end square
        //and drawing the stroke
        ctx.beginPath();
        ctx.moveTo(fromx, fromy);
        ctx.lineTo(tox, toy);
        ctx.lineWidth = arrowWidth;
        ctx.stroke();

        //starting a new path from the head of the arrow to one of the sides of
        //the point
        ctx.beginPath();
        ctx.moveTo(tox, toy);
        ctx.lineTo(tox - headlen * Math.cos(angle - Math.PI / 7),
            toy - headlen * Math.sin(angle - Math.PI / 7));

        //path from the side point of the arrow, to the other side point
        ctx.lineTo(tox - headlen * Math.cos(angle + Math.PI / 7),
            toy - headlen * Math.sin(angle + Math.PI / 7));

        //path from the side point back to the tip of the arrow, and then
        //again to the opposite side point
        ctx.lineTo(tox, toy);
        ctx.lineTo(tox - headlen * Math.cos(angle - Math.PI / 7),
            toy - headlen * Math.sin(angle - Math.PI / 7));

        //draws the paths created above
        ctx.stroke();
        ctx.restore();
    }
}
