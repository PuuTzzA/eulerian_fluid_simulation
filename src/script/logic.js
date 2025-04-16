const canvas = document.getElementById("canvas");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
const ctx = canvas.getContext("2d");

const FLUID = 1;
const SOLID = 2;
const EMPTY = 3;

window.addEventListener("resize", e => {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
})

class GraphicsSettings {
    constructor() {
        this.lineWidth = 1;
        this.lineColour = "transparent";
        this.particleSize = 3;
        this.maxSpeed = 10;
        this.maxVelocitySize = 1;
        this.arrowWidth = 1;
        this.arrowHeadSize = 9;
        this.arrowColor = "red";
    }

    getLineWidth() {
        return this.lineWidth;
    }

    getLineColour() {
        return this.lineColour;
    }

    getHatchingSettings(label, gridSize) {
        const solidSettings = {
            colour: "rgb(112, 101, 84)",
            lineWidth: 3,
            cross: true,
            direction: "right",
            delta: gridSize / 6
        };
        const fluidSettings = {
            colour: "rgb(28, 164, 255)",
            lineWidth: 1,
            cross: false,
            direction: "left",
            delta: gridSize / 4
        };
        const emptySettings = {
            colour: "white",
            lineWidth: 0.5,
            cross: false,
            direction: "right",
            delta: gridSize / 2
        };

        switch (label) {
            case FLUID:
                return fluidSettings;
            case SOLID:
                return solidSettings;
            case EMPTY:
                return emptySettings;
            default:
                return "rgb(255, 0, 255)";
        }
    }

    getHatchingDelta() {

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

class FluidQuantity {
    constructor(w, h, offsetX, offsetY, cellSize) {
        this.w = w;
        this.h = h;
        this.offsetX = offsetX;
        this.offsetY = offsetY;
        this.cellSize = cellSize;

        this.src = new Float64Array(w * h);
        this.dst = new Float64Array(w * h);
    }

    flip() {
        const tmp = this.src;
        this.src = this.dst;
        this.dst = tmp;
    }

    at(x, y) {
        return this.src[y * this.w + x];
    }

    id(x, y) {
        return y * this.w + x;
    }

    lerp(x, y, t) {
        return x * (1 - t) + y * t;
    }

    cerp(a, b, c, d, x) {
        const xsq = x * x;
        const xcu = xsq * x;

        const minV = Math.min(a, Math.min(b, Math.min(c, d)));
        const maxV = Math.max(a, Math.max(b, Math.max(c, d)));

        const t =
            a * (-(1 / 3) * x + 0.5 * xsq - (1 / 6) * xcu) +
            b * (1 - xsq + 0.5 * (xcu - x)) +
            c * (x + 0.5 * (xsq - xcu)) +
            d * ((1 / 6) * (xcu - x));

        return Math.min(Math.max(t, minV), maxV);
    }

    bilinearInterpolation(x, y) {
        x = Math.min(Math.max(x - this.offsetX, 0), this.w - 1.00001);
        y = Math.min(Math.max(y - this.offsetY, 0), this.h - 1.00001);

        const idX = Math.floor(x);
        const idY = Math.floor(y);

        x -= idX;
        y -= idY;

        const x00 = this.at(idX + 0, idY + 0);
        const x10 = this.at(idX + 1, idY + 0);
        const x01 = this.at(idX + 0, idY + 1);
        const x11 = this.at(idX + 1, idY + 1);

        return this.lerp(this.lerp(x00, x10, x), this.lerp(x01, x11, x), y);
    }

    bicubicInterpolation(x, y) {
        x = Math.min(Math.max(x - this.offsetX, 0), this.w - 1.00001);
        y = Math.min(Math.max(y - this.offsetY, 0), this.h - 1.00001);

        const ix = Math.floor(x);
        const iy = Math.floor(y);

        x -= ix;
        y -= iy;

        const x0 = Math.max(ix - 1, 0);
        const x1 = ix;
        const x2 = ix + 1;
        const x3 = Math.min(ix + 2, this.w - 1);

        const y0 = Math.max(iy - 1, 0);
        const y1 = iy;
        const y2 = iy + 1;
        const y3 = Math.min(iy + 2, this.h - 1);

        const q0 = this.cerp(this.at(x0, y0), this.at(x1, y0), this.at(x2, y0), this.at(x3, y0), x);
        const q1 = this.cerp(this.at(x0, y1), this.at(x1, y1), this.at(x2, y1), this.at(x3, y1), x);
        const q2 = this.cerp(this.at(x0, y2), this.at(x1, y2), this.at(x2, y2), this.at(x3, y2), x);
        const q3 = this.cerp(this.at(x0, y3), this.at(x1, y3), this.at(x2, y3), this.at(x3, y3), x);

        return this.cerp(q0, q1, q2, q3, y);
    }

    explicitEuler(x, y, dt, velocityX, velocityY) {
        const uVel = velocityX.bilinearInterpolation(x, y) / this.cellSize;
        const vVel = velocityY.bilinearInterpolation(x, y) / this.cellSize;

        x -= uVel * dt;
        y -= vVel * dt;

        return [x, y];
    }

    rungeKutta4(x, y, dt, velocityX, velocityY) {
        const uk1 = velocityX.bilinearInterpolation(x, y) / this.cellSize;
        const vk1 = velocityY.bilinearInterpolation(x, y) / this.cellSize;
        const xk1 = x - 0.5 * dt * uk1;
        const yk1 = y - 0.5 * dt * vk1;

        const uk2 = velocityX.bilinearInterpolation(xk1, yk1) / this.cellSize;
        const vk2 = velocityY.bilinearInterpolation(xk1, yk1) / this.cellSize;
        const xk2 = x - 0.5 * dt * uk2;
        const yk2 = y - 0.5 * dt * vk2;

        const uk3 = velocityX.bilinearInterpolation(xk2, yk2) / this.cellSize;
        const vk3 = velocityY.bilinearInterpolation(xk2, yk2) / this.cellSize;
        const xk3 = x - dt * uk3;
        const yk3 = y - dt * vk3;

        const uk4 = velocityX.bilinearInterpolation(xk3, yk3) / this.cellSize;
        const vk4 = velocityY.bilinearInterpolation(xk3, yk3) / this.cellSize;

        x -= dt / 6 * (uk1 + 2 * uk2 + 2 * uk3 + uk4);
        y -= dt / 6 * (vk1 + 2 * vk2 + 2 * vk3 + vk4);

        return [x, y];
    }

    advect(dt, u, v) {
        for (let iy = 0, idx = 0; iy < this.h; iy++) {
            for (let ix = 0; ix < this.w; ix++, idx++) {
                let x = ix + this.offsetX;
                let y = iy + this.offsetY;

                [x, y] = this.rungeKutta4(x, y, dt, u, v);

                this.dst[idx] = this.bicubicInterpolation(x, y);
            }
        }
    }

    static cubicPulse(x) {
        x = Math.min(Math.abs(x), 1);
        return 1 - x * x * (3 - 2 * x);
    }

    static length(x, y) {
        return Math.sqrt(x * x + y * y);
    }

    addInflow(x0, y0, x1, y1, value) {
        const ix0 = Math.floor(x0 / this.cellSize - this.offsetX);
        const iy0 = Math.floor(y0 / this.cellSize - this.offsetY);
        const ix1 = Math.floor(x1 / this.cellSize - this.offsetX);
        const iy1 = Math.floor(y1 / this.cellSize - this.offsetY);

        for (let y = Math.max(iy0, 0); y < Math.min(iy1, this.h); y++) {
            for (let x = Math.max(ix0, 0); x < Math.min(ix1, this.w); x++) {

                const l = FluidQuantity.length(
                    (2 * (x + 0.5) * this.cellSize - (x0 + x1)) / (x1 - x0),
                    (2 * (y + 0.5) * this.cellSize - (y0 + y1)) / (y1 - y0)
                )

                const vi = FluidQuantity.cubicPulse(l) * value;

                if (Math.abs(this.src[y * this.w + x]) < Math.abs(vi)) {
                    this.src[y * this.w + x] = vi;
                }
            }
        }
    }
}

class Fluid {
    constructor(w, h, density) {
        this.gs = new GraphicsSettings();

        this.w = w;
        this.h = h;
        this.density = density;
        this.cellSize = 1 / Math.min(w, h);

        this.d = new FluidQuantity(w, h, 0.5, 0.5, this.cellSize);
        this.u = new FluidQuantity(w + 1, h, 0.0, 0.5, this.cellSize);
        this.v = new FluidQuantity(w, h + 1, 0.5, 0.0, this.cellSize);

        this.rhs = new Float64Array(w * h);
        this.pressure = new Float64Array(w * h);


        const shouldX = window.innerWidth / w;
        const shouldY = window.innerHeight / h;
        this.gridPixelSize = Math.min(shouldX, shouldY);
    }

    update(dt) {
        this.buildRhs(dt);
        this.project(600, dt);
        this.applyPressure(dt);

        this.d.advect(dt, this.u, this.v);
        this.u.advect(dt, this.u, this.v);
        this.v.advect(dt, this.u, this.v);

        this.d.flip();
        this.u.flip();
        this.v.flip();
    }

    addInflow(x, y, w, h, d, u, v) {
        this.d.addInflow(x, y, x + w, y + h, d);
        this.u.addInflow(x, y, x + w, y + h, u);
        this.v.addInflow(x, y, x + w, y + h, v);
    }

    maxTimestep() {
        let maxVel = 0;
        for (let y = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++) {

                const u = this.u.bilinearInterpolation(x + 0.5, y + 0.5);
                const v = this.v.bilinearInterpolation(x + 0.5, y + 0.5);

                const vel = Math.sqrt(u * u + v * v);
                maxVel = Math.max(maxVel, vel);
            }
        }

        // Fluid should not move more than two grid cells per iteration
        const maxTimestep = 0.5 * this.cellSize / maxVel;
        return Math.min(maxTimestep, 1);
    }

    draw() {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        ctx.lineWidth = this.gs.getLineWidth();

        const offsetTop = window.innerHeight - this.gridPixelSize * this.h;

        for (let y = 0; y < this.d.h; y++) {
            for (let x = 0; x < this.d.w; x++) {

                ctx.strokeStyle = this.gs.getLineColour();
                const d = Math.floor(this.d.at(x, y) * 100);
                ctx.fillStyle = this.gs.getCellColour(d, 0, 0);

                ctx.strokeRect(x * this.gridPixelSize, y * this.gridPixelSize + offsetTop, this.gridPixelSize, this.gridPixelSize);
                ctx.fillRect(x * this.gridPixelSize, y * this.gridPixelSize + offsetTop, this.gridPixelSize, this.gridPixelSize);

                // Fluid.hatchRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize, this.gs.getHatchingSettings(this.labels[i], this.gridSize));
            }
        }

        return;

        for (let y = 0; y < this.u.h; y++) {
            for (let x = 0; x < this.u.w; x++) {

                if (Math.abs(this.u.at(x, y)) < 0.01) continue;

                let fromX = x * this.gridPixelSize;
                let fromY = (y + 0.5) * this.gridPixelSize + offsetTop;
                let toX = (x + this.gs.getArrowLength(this.u.at(x, y))) * this.gridPixelSize;
                let toY = fromY;

                let col = this.gs.getArrowColor();
                Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
            }
        }

        for (let y = 0; y < this.v.h; y++) {
            for (let x = 0; x < this.v.w; x++) {

                if (Math.abs(this.v.at(x, y)) < 0.01) continue;

                let fromX = (x + 0.5) * this.gridPixelSize;
                let fromY = y * this.gridPixelSize + offsetTop;
                let toX = fromX;
                let toY = (y + this.gs.getArrowLength(this.v.at(x, y))) * this.gridPixelSize + offsetTop;

                let col = this.gs.getArrowColor();

                Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
            }
        }
    }

    buildRhs() {
        // builds the pressure right hand side, meaning the negative divergence
        const scale = 1 / this.cellSize;

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                const du = this.u.at(x + 1, y) - this.u.at(x, y);
                const dv = this.v.at(x, y + 1) - this.v.at(x, y);
                this.rhs[idx] = -scale * (du + dv);
            }
        }
    }

    project(limit, dt) {
        // gauss seidel iteration
        const scale = dt / (this.density * this.cellSize * this.cellSize);

        let maxDelta;

        for (let iter = 0; iter < limit; iter++) {
            maxDelta = 0;

            for (let y = 0, idx = 0; y < this.h; y++) {
                for (let x = 0; x < this.w; x++, idx++) {
                    let diag = 0;
                    let offDiag = 0;

                    if (x > 0) {
                        diag += scale;
                        offDiag += scale * this.pressure[idx - 1];
                    }
                    if (y > 0) {
                        diag += scale;
                        offDiag += scale * this.pressure[idx - this.w];
                    }
                    if (x < this.w - 1) {
                        diag += scale;
                        offDiag += scale * this.pressure[idx + 1];
                    }
                    if (y < this.h - 1) {
                        diag += scale;
                        offDiag += scale * this.pressure[idx + this.w];
                    }

                    const newPressure = (this.rhs[idx] + offDiag) / diag;
                    maxDelta = Math.max(maxDelta, Math.abs(newPressure - this.pressure[idx]));
                    this.pressure[idx] = newPressure;
                }
            }

            if (maxDelta < 1e-5) {
                console.log(`Exiting solver after ${iter} iterations`)
                return;
            }
        }

        console.log("EXCEEDED MAXIMUM ITERATIONS");
    }

    applyPressure(dt) {
        const scale = dt / (this.density * this.cellSize);

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {

                this.u.src[this.u.id(x, y)] -= scale * this.pressure[idx];
                this.u.src[this.u.id(x + 1, y)] += scale * this.pressure[idx];

                this.v.src[this.v.id(x, y)] -= scale * this.pressure[idx];
                this.v.src[this.v.id(x, y + 1)] += scale * this.pressure[idx];
            }
        }

        for (let y = 0; y < this.h; y++) {
            this.u.src[y * (this.w + 1)] = 0;
            this.u.src[y * (this.w + 1) + this.w] = 0;
        }
        for (let x = 0; x < this.w; x++) {
            this.v.src[x] = 0;
            this.v.src[this.h * this.w + x] = 0;
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

    static hatchRect(x1, y1, dx, dy, settings) {
        const cross = settings.cross;
        const direction = settings.direction;
        const color = settings.colour;
        const delta = settings.delta;
        const lw = settings.lineWidth;

        ctx.save();
        ctx.beginPath();
        ctx.rect(x1, y1, dx, dy);
        ctx.clip();

        const majorAxe = Math.max(dx, dy);
        ctx.strokeStyle = color;
        ctx.lineWidth = lw;

        const drawLines = (dir) => {
            for (let n = -majorAxe; n < majorAxe; n += delta) {
                ctx.beginPath();
                if (dir === 'right') {
                    // bottom-left to top-right (like your original)
                    ctx.moveTo(n + x1, y1);
                    ctx.lineTo(n + x1 + dy, y1 + dy);
                } else if (dir === 'left') {
                    // top-left to bottom-right
                    ctx.moveTo(n + x1, y1 + dy);
                    ctx.lineTo(n + x1 + dy, y1);
                }
                ctx.stroke();
            }
        };

        // Draw primary direction
        drawLines(direction);

        // Optionally add crosshatching
        if (cross) {
            drawLines(direction === 'right' ? 'left' : 'right');
        }

        ctx.restore();
    }
}
