const canvas = document.getElementById("canvas");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
const ctx = canvas.getContext("2d");

const CELL_FLUID = 0;
const CELL_SOLID = 1;
const CELL_EMPTY = 2;

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

/**
 * -1, 0 or 1
 */
function signum(x) {
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

/**
 * -1 or 1
 */
function nonZeroSgn(x) {
    return x < 0 ? -1 : 1;
}

function cubicPulse(x) {
    x = Math.min(Math.abs(x), 1);
    return 1 - x * x * (3 - 2 * x);
}

function length(x, y) {
    return Math.sqrt(x * x + y * y);
}

function scaledAdd(dst, a, b, s, cell) {
    // 'dst' = 'a' + s * 'b'
    for (let i = 0; i < dst.length; i++) {
        if (cell[i] == CELL_FLUID) {
            dst[i] = a[i] + s * b[i];
        }
    }
}

function infinityNorm(a, cell) {
    // returns the infinity norm of vector 'a'
    let max = 0;
    for (let i = 0; i < a.length; i++) {
        if (cell[i] == CELL_FLUID) {
            max = Math.max(max, Math.abs(a[i]));
        }
    }
    return max;
}

function dotProduct(a, b, cell) {
    let sum = 0;
    for (let i = 0; i < a.length; i++) {
        if (cell[i] == CELL_FLUID) {
            sum += a[i] * b[i];
        }
    }
    return sum;
}

/**
 * rotates point (x, y) by angle phi
 */
function rotate(x, y, phi) {
    const tmpX = x, tmpY = y;
    x = Math.cos(phi) * tmpX + Math.sin(phi) * tmpY;
    y = -Math.sin(phi) * tmpX + Math.cos(phi) * tmpY;
    return [x, y];
}

/**
 * For three corners in a 1x1 square, with 'in' being inside the fluid, and adjacent to 
 * to 'out1' and 'out2'. All three values are distances to a surface.
 * Returns fraction of the full square occupied by the surface in the specific triangle.
 */
function triangleOccupancy(out1, in1, out2) {
    return 0.5 * in1 * in1 / ((out1 - in1) * (out2 - in1));
}

/**
 * For four courners in a 1x1 square (two inside, two outside)
 * Returns fraction of the square occupied by the surface.
 */
function trapezoidOccupancy(out1, out2, in1, in2) {
    return 0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2));
}

/**
 * given the distance of four corner points of a 1x1 square, returns the 
 * area of the part of the square occupied by the surface.
 */
function occupancy(d11, d12, d21, d22) {
    const ds = [d11, d12, d22, d21];
    let mask = 0;
    for (let i = 3; i >= 0; i--) {
        mask = (mask << 1) | (ds[i] < 0 ? 1 : 0);
    }
    // mask now is D21 D22 D12 D11 whre DXX is 0 if inside, 1 if outside
    switch (mask) {
        case 0b0000: // all outside
            return 0;
        // one inside
        case 0b0001:
            return triangleOccupancy(d21, d11, d12);
        case 0b0010:
            return triangleOccupancy(d11, d12, d22);
        case 0b0100:
            return triangleOccupancy(d12, d22, d21);
        case 0b1000:
            return triangleOccupancy(d22, d21, d11);
        // one outside
        case 0b1110:
            return 1 - triangleOccupancy(-d21, -d11, -d12);
        case 0b1101:
            return 1 - triangleOccupancy(-d11, -d12, -d22);
        case 0b1011:
            return 1 - triangleOccupancy(-d12, -d22, -d21);
        case 0b0111:
            return 1 - triangleOccupancy(-d22, -d21, -d11);
        // two adjacent inside
        case 0b0011:
            return trapezoidOccupancy(d21, d22, d11, d12);
        case 0b0110:
            return trapezoidOccupancy(d11, d21, d12, d22);
        case 0b1100:
            return trapezoidOccupancy(d11, d12, d21, d22);
        case 0b1001:
            return trapezoidOccupancy(d12, d22, d11, d21);
        // two opposed inside
        case 0b0101:
            return triangleOccupancy(d21, d11, d12) + triangleOccupancy(d12, d22, d21);
        case 0b1010:
            return triangleOccupancy(d11, d12, d22) + triangleOccupancy(d12, d22, d21);
        // all inside
        case 0b1111:
            return 1;
    }
}

/**
 * abstract class for solid body
 */
class SolidBody {
    constructor(posX, posY, scaleX, scaleY, theta, velX, velY, velTheta) {
        if (this.constructor == SolidBody) {
            throw new Error("Abstract classes can't be instantiated.");
        }

        this.posX = posX;
        this.posY = posY;
        this.scaleX = scaleX;
        this.scaleY = scaleY;
        this.theta = theta;
        this.velX = velX;
        this.velY = velY;
        this.velTheta = velTheta;
    }

    /**
     * transforms point(x, y) from the global to the local coordinate system
     */
    globalToLocal(x, y) {
        x -= this.posX;
        y -= this.posY;
        [x, y] = rotate(x, y, -this.theta);
        x /= this.scaleX;
        y /= this.scaleY;
        return [x, y];
    }

    /**
     * transforms point(x, y) from the local to the global coordinate system
     */
    localToGlobal(x, y) {
        x *= this.scaleX;
        y *= this.scaleY;
        [x, y] = rotate(x, y, this.theta);
        x += this.posX;
        y += this.posY;
        return [x, y];
    }

    /**
     * returns the signed distance from (x, y) to the closest point on the solids surface
     */
    distance(x, y) {
        throw new Error("Method must be implemented.");
    }

    /**
     * returns the closest point on the surface of the solid 
     */
    closestSurfacePoint(x, y) {
        throw new Error("Method must be implemented.");
    }

    /**
     * returns the gradient at point (x, y)
     */
    distanceNormal(x, y) {
        throw new Error("Method must be implemented.");
    }

    velocityX(x, y) {
        x -= this.posX;
        y -= this.posY;
        return this.velX - y * this.velTheta;
    }

    velocityY(x, y) {
        x -= this.posX;
        y -= this.posY;
        return this.velY + x * this.velTheta;
    }

    velocity(x, y) {
        return [this.velocityX(x, y), this.velocityY(x, y)];
    }

    update(dt) {
        this.posX += this.velX * dt;
        this.posY += this.velY * dt;
        this.theta += this.velTheta * dt;
    }
}

class SolidBox extends SolidBody {
    constructor(posX, posY, scaleX, scaleY, theta, velX, velY, velT) {
        super(posX, posY, scaleX, scaleY, theta, velX, velY, velT);
    }

    distance(x, y) {
        x -= this.posX;
        y -= this.posY;
        [x, y] = rotate(x, y, -this.theta);
        const dx = Math.abs(x) - this.scaleX * 0.5;
        const dy = Math.abs(y) - this.scaleY * 0.5;

        if (dx >= 0 || dy >= 0) {
            return length(Math.max(dx, 0), Math.max(dy, 0));
        }
        return Math.max(dx, dy);
    }

    closestSurfacePoint(x, y) {
        x -= this.posX;
        y -= this.posY;
        [x, y] = rotate(x, y, -this.theta);
        const dx = Math.abs(x) - this.scaleX * 0.5;
        const dy = Math.abs(y) - this.scaleY * 0.5;

        if (dx > dy) {
            x = nonZeroSgn(x) * 0.5 * this.scaleX;
        } else {
            y = nonZeroSgn(y) * 0.5 * this.scaleY;
        }

        if (dx >= 0 && dy >= 0) {
            x = nonZeroSgn(x) * 0.5 * this.scaleX;
            y = nonZeroSgn(y) * 0.5 * this.scaleY;
        }

        [x, y] = rotate(x, y, this.theta);
        x += this.posX;
        y += this.posY;
        return [x, y];
    }

    distanceNormal(x, y) {
        x -= this.posX;
        y -= this.posY;
        [x, y] = rotate(x, y, -this.theta);

        let nx, ny;
        if (Math.abs(x) - this.scaleX * 0.5 > Math.abs(y) - this.scaleY * 0.5) {
            nx = nonZeroSgn(x);
            ny = 0;
        } else {
            nx = 0;
            ny = nonZeroSgn(y);
        }

        [nx, ny] = rotate(nx, ny, this.theta);
        return [nx, ny];
    }
}

class SolidSphere extends SolidBody {
    constructor(posX, posY, scale, theta, velX, velY, velT) {
        super(posX, posY, scale, scale, theta, velX, velY, velT);
    }

    distance(x, y) {
        return length(x - this.posX, y - this.posY) - this.scaleX * 0.5;
    }

    closestSurfacePoint(x, y) {
        [x, y] = this.globalToLocal(x, y);

        const r = length(x, y);

        if (r < 0.0001) {
            x = 0.5;
            y = 0.0;
        } else {
            x /= 2 * r;
            y /= 2 * r
        }

        [x, y] = this.localToGlobal(x, y);
        return [x, y];
    }

    distanceNormal(x, y) {
        x -= this.posX;
        y -= this.posY;
        const r = length(x, y);

        let nx, ny;
        if (r < 0.0001) {
            nx = 1;
            ny = 0;
        } else {
            nx = x / r;
            ny = y / r;
        }
        return [nx, ny];
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

        this.phi = new Float64Array((w + 1) * (h + 1)); // distance field, samples are offset by (-0.5, -0.5) and the grid is larger in each dimension, so that every src has 4 phi values 
        this.volume = new Float64Array(w * h); // fractional cell volume occupied by the fluid

        this.normalX = new Float64Array(w * h);
        this.normalY = new Float64Array(w * h);

        this.cell = new Uint8Array(w * h); // Fluid or solid cell (only solid if cell completely covered by bodies)
        this.body = new Uint8Array(w * h); // specifies the index of the solid body closest to a grid cell
        this.mask = new Uint8Array(w * h); // auxiliary array for extrapolation

        this.cell.fill(CELL_FLUID);
        this.volume.fill(1);
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

    volumeAt(x, y) {
        return this.volume[y * this.w + x];
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

    /**
     * if (x, y) is inside a solid, project it back out to the closest surface point
     */
    backProject(x, y, bodies) {
        const rx = Math.min(Math.max(Math.floor(x - this.offsetX), 0), this.w - 1);
        const ry = Math.min(Math.max(Math.floor(y - this.offsetY), 0), this.h - 1);

        if (this.cell[ry * this.w + rx] != CELL_FLUID) {
            x = (x - this.offsetX) * this.cellSize;
            y = (y - this.offsetY) * this.cellSize;

            [x, y] = bodies[this.body[ry * this.w + rx]].closestSurfacePoint(x, y);
            x = x / this.cellSize + this.offsetX;
            y = y / this.cellSize + this.offsetY;
        }
        return [x, y];
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

    advect(dt, u, v, bodies) {
        for (let iy = 0, idx = 0; iy < this.h; iy++) {
            for (let ix = 0; ix < this.w; ix++, idx++) {
                let x = ix + this.offsetX;
                let y = iy + this.offsetY;

                [x, y] = this.rungeKutta4(x, y, dt, u, v);
                [x, y] = this.backProject(x, y, bodies);

                this.dst[idx] = this.bicubicInterpolation(x, y);
            }
        }
    }

    addInflow(x0, y0, x1, y1, value) {
        const ix0 = Math.floor(x0 / this.cellSize - this.offsetX);
        const iy0 = Math.floor(y0 / this.cellSize - this.offsetY);
        const ix1 = Math.floor(x1 / this.cellSize - this.offsetX);
        const iy1 = Math.floor(y1 / this.cellSize - this.offsetY);

        for (let y = Math.max(iy0, 0); y < Math.min(iy1, this.h); y++) {
            for (let x = Math.max(ix0, 0); x < Math.min(ix1, this.w); x++) {

                const l = length(
                    (2 * (x + 0.5) * this.cellSize - (x0 + x1)) / (x1 - x0),
                    (2 * (y + 0.5) * this.cellSize - (y0 + y1)) / (y1 - y0)
                )

                const vi = cubicPulse(l) * value;

                if (Math.abs(this.src[y * this.w + x]) < Math.abs(vi)) {
                    this.src[y * this.w + x] = vi;
                }
            }
        }
    }

    /**
     * fills all solid related fields (phi, volume, cell, body and normalX/Y)
     */
    fillSolidFields(bodies) {
        if (bodies.length == 0) return;

        // compute distance field first
        for (let iy = 0, idx = 0; iy <= this.h; iy++) {
            for (let ix = 0; ix <= this.w; ix++, idx++) {
                const x = (ix + this.offsetX - 0.5) * this.cellSize;
                const y = (iy + this.offsetY - 0.5) * this.cellSize;

                this.phi[idx] = bodies[0].distance(x, y);
                for (let i = 1; i < bodies.length; i++) {
                    this.phi[idx] = Math.min(this.phi[idx], bodies[i].distance(x, y));
                }
            }
        }

        for (let iy = 0, idx = 0; iy < this.h; iy++) {
            for (let ix = 0; ix < this.w; ix++, idx++) {
                const x = (ix + this.offsetX) * this.cellSize;
                const y = (iy + this.offsetY) * this.cellSize;

                this.body[idx] = 0;
                let d = bodies[0].distance(x, y);

                for (let i = 1; i < bodies.length; i++) {
                    const id = bodies[i].distance(x, y);
                    if (id < d) {
                        this.body[idx] = i;
                        d = id;
                    }
                }

                // compute cell volume from the four adjacent distance values
                const idxp = iy * (this.w + 1) + ix;
                this.volume[idx] = 1 - occupancy(
                    this.phi[idxp], this.phi[idxp + 1],
                    this.phi[idxp + this.w + 1], this.phi[idxp + this.w + 2]
                );

                if (this.volume[idx] < 0.01) {
                    this.volume[idx] = 0;
                }

                let nx, ny;
                [nx, ny] = bodies[this.body[idx]].distanceNormal(x, y);
                this.normalX[idx] = nx;
                this.normalY[idx] = ny;

                this.cell[idx] = this.volume[idx] == 0 ? CELL_SOLID : CELL_FLUID;
            }
        }
    }

    /**
     * If an entry is 0, it means both neighbours are available and the cell is ready for the PDE solve. 
     * If it is 1, the cell waits for neighbour in x direction, 2 for y direction and 3 for both.
     */
    fillSolidMask() {
        for (let y = 1; y < this.h - 1; y++) {
            for (let x = 1; x < this.w - 1; x++) {
                const idx = y * this.w + x;

                if (this.cell[idx] == CELL_FLUID) continue;

                let nx = this.normalX[idx];
                let ny = this.normalY[idx];

                this.mask[idx] = 0;
                if (nx != 0 && this.cell[idx + signum(nx)] != CELL_FLUID) {
                    this.mask[idx] |= 1; // neighbour in x direction is blocked
                }
                if (ny != 0 && this.cell[idx + this.w * signum(ny)] != CELL_FLUID) {
                    this.mask[idx] |= 2; // neighbour in y direction is blocked
                }
            }
        }
    }

    /**
     * Solve for value at index idx.
     * The value is computed such that the directional derivative along the disstance field normal is 0.
     */
    extrapolateNormal(idx) {
        const nx = this.normalX[idx];
        const ny = this.normalY[idx];

        const srcX = this.src[idx + signum(nx)];
        const srcY = this.src[idx + this.w * signum(ny)];

        return (Math.abs(nx) * srcX + Math.abs(ny) * srcY) / (Math.abs(nx) + Math.abs(ny));
    }

    /**
     * Given that a neighbour (specified by mask [1 = x, 2 = y]) has been solved, update the mask 
     * and if the cell can now be computed, add it to the queue of ready cells.
     */
    freeNeighbour(idx, border, mask) {
        this.mask[idx] &= ~mask; // remove the mask bit
        if (this.cell[idx] != CELL_FLUID && this.mask[idx] == 0) {
            border.push(idx);
        }
    }

    /**
     * Extrapolate fluid quantities into solids, where these quantities would normally be undefined.
     * Set quantities to the closest value on the solid-fluid boundary.
     * Or such that the gradient of the quantity is 0 along the gradient of the distance field. 
     */
    extrapolate() {
        this.fillSolidMask();

        let border = [];

        for (let y = 1; y < this.h - 1; y++) {
            for (let x = 1; x < this.w - 1; x++) {
                const idx = y * this.w + x;

                if (this.cell[idx] != CELL_FLUID && this.mask[idx] == 0) {
                    border.push(idx);
                }
            }
        }

        while (border.length != 0) {
            const idx = border.pop();

            // solve for the value in cell
            this.src[idx] = this.extrapolateNormal(idx);

            // Notify adjacent cells
            if (this.normalX[idx - 1] > 0) {
                this.freeNeighbour(idx - 1, border, 1);
            }
            if (this.normalX[idx + 1] < 0) {
                this.freeNeighbour(idx + 1, border, 1);
            }
            if (this.normalY[idx - this.w] > 0) {
                this.freeNeighbour(idx - this.w, border, 2);
            }
            if (this.normalY[idx + this.w] < 0) {
                this.freeNeighbour(idx + this.w, border, 2);
            }
        }
    }
}

class Fluid {
    constructor(w, h, densityAir, densitySoot, diffusion, bodies) {
        this.gs = new GraphicsSettings();
        this.mouseX = 0;
        this.mouseY = 0;

        this.w = w;
        this.h = h;
        this.densityAir = densityAir;
        this.densitySoot = densitySoot;
        this.diffusion = diffusion;
        this.cellSize = 1 / Math.min(w, h);

        this.tAmb = 294; // ambient temperature in Kelvin
        this.gravityX = 0;
        this.gravityY = 9.81;

        this.d = new FluidQuantity(w, h, 0.5, 0.5, this.cellSize);
        this.t = new FluidQuantity(w, h, 0.5, 0.5, this.cellSize);
        this.u = new FluidQuantity(w + 1, h, 0.0, 0.5, this.cellSize);
        this.v = new FluidQuantity(w, h + 1, 0.5, 0.0, this.cellSize);

        this.t.src.fill(this.tAmb);

        this.rhs = new Float64Array(w * h);
        this.pressure = new Float64Array(w * h);

        this.z = new Float64Array(w * h);
        this.s = new Float64Array(w * h);
        this.aDiag = new Float64Array(w * h);
        this.aPlusX = new Float64Array(w * h);
        this.aPlusY = new Float64Array(w * h);
        this.precon = new Float64Array(w * h);

        this.uDensity = new Float64Array((w + 1) * h);
        this.vDensity = new Float64Array(w * (h + 1));

        this.bodies = bodies;

        const shouldX = window.innerWidth / w;
        const shouldY = window.innerHeight / h;
        this.gridPixelSize = Math.min(shouldX, shouldY);
    }

    update(dt) {
        this.d.fillSolidFields(this.bodies);
        this.t.fillSolidFields(this.bodies);
        this.u.fillSolidFields(this.bodies);
        this.v.fillSolidFields(this.bodies);

        // Temperature diffusion
        this.rhs.set(this.t.src);
        this.buildHeatDiffusionMatrix(dt);
        this.buildPreconditioner();
        this.projectConjugateGradient(2000);
        this.t.src.set(this.pressure);
        this.t.extrapolate();

        this.addBuoyancy(dt);
        this.setBoundaryCondition();

        // Incompressility 
        this.buildRhs(dt);
        this.computeDensities();
        this.buildPressureMatrix(dt);
        this.buildPreconditioner();
        this.projectConjugateGradient(2000);
        this.applyPressure(dt);

        this.d.extrapolate();
        this.u.extrapolate();
        this.v.extrapolate();

        this.setBoundaryCondition();

        this.d.advect(dt, this.u, this.v, this.bodies);
        this.t.advect(dt, this.u, this.v, this.bodies);
        this.u.advect(dt, this.u, this.v, this.bodies);
        this.v.advect(dt, this.u, this.v, this.bodies);

        this.d.flip();
        this.t.flip();
        this.u.flip();
        this.v.flip();
    }

    addInflow(x, y, w, h, d, t, u, v) {
        this.d.addInflow(x, y, x + w, y + h, d);
        this.t.addInflow(x, y, x + w, y + h, t);
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

    /**
     * Builds the pressure right hand side, meaning the negative divergence.
     * "Blend" between solid and fluid velocity based on the cell volume.
     */
    buildRhs() {
        const cell = this.d.cell;
        const body = this.d.body;

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] == CELL_FLUID) {

                    let e = 0;
                    if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                        e += this.u.volumeAt(x, y) * this.u.at(x, y);
                        e += (1 - this.u.volumeAt(x, y)) * this.bodies[this.u.body[this.u.id(x, y)]].velocityX((x - 0.5) * this.cellSize, y * this.cellSize);
                    }
                    if (y > 0 && cell[idx - this.w] == CELL_FLUID) {
                        e += this.v.volumeAt(x, y) * this.v.at(x, y);
                        e += (1 - this.v.volumeAt(x, y)) * this.bodies[this.v.body[this.v.id(x, y)]].velocityY(x * this.cellSize, (y - 0.5) * this.cellSize);
                    }
                    if (x < this.w - 1 && cell[idx + 1] == CELL_FLUID) {
                        e -= this.u.volumeAt(x + 1, y) * this.u.at(x + 1, y);
                        e -= (1 - this.u.volumeAt(x + 1, y)) * this.bodies[this.u.body[this.u.id(x + 1, y)]].velocityX((x + 0.5) * this.cellSize, y * this.cellSize);
                    }
                    if (y < this.h - 1 && cell[idx + this.w] == CELL_FLUID) {
                        e -= this.v.volumeAt(x, y + 1) * this.v.at(x, y + 1);
                        e -= (1 - this.v.volumeAt(x, y + 1)) * this.bodies[this.v.body[this.v.id(x, y + 1)]].velocityY(x * this.cellSize, (y + 0.5) * this.cellSize);
                    }

                    this.rhs[idx] = e;
                } else {
                    this.rhs[idx] = 0;
                }
            }
        }
    }

    /**
     * Computes densities as a function of temperature and smoke concentration.
     */
    computeDensities() {
        const alpha = (this.densitySoot - this.densityAir) / this.densityAir;

        this.uDensity.fill(0);
        this.vDensity.fill(0);

        for (let y = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++) {
                let density = this.densityAir * this.tAmb / this.t.at(x, y) * (1 + alpha * this.d.at(x, y));
                density = Math.max(density, 0.05 * this.densityAir);

                this.uDensity[this.u.id(x, y)] += 0.5 * density;
                this.vDensity[this.v.id(x, y)] += 0.5 * density;
                this.uDensity[this.u.id(x + 1, y)] += 0.5 * density;
                this.vDensity[this.v.id(x, y + 1)] += 0.5 * density;
            }
        }
    }

    /**
     * Builds the pressure matrix (A)
     */
    buildPressureMatrix(dt) {
        const scale = dt / (this.cellSize);
        const cell = this.d.cell;

        this.aDiag.fill(0);
        this.aPlusX.fill(0);
        this.aPlusY.fill(0);

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] != CELL_FLUID) {
                    continue;
                }
                if (x < this.w - 1 && cell[idx + 1] == CELL_FLUID) {
                    const factor = scale * this.u.volumeAt(x + 1, y) / this.uDensity[this.u.id(x + 1, y)];
                    this.aDiag[idx] += factor;
                    this.aDiag[idx + 1] += factor;
                    this.aPlusX[idx] = -factor;
                }
                if (y < this.h - 1 && cell[idx + this.w] == CELL_FLUID) {
                    const factor = scale * this.v.volumeAt(x, y + 1) / this.vDensity[this.v.id(x, y + 1)];
                    this.aDiag[idx] += factor;
                    this.aDiag[idx + this.w] += factor;
                    this.aPlusY[idx] = -factor;
                }
            }
        }
    }

    /**
     * Builds the matrix used to compute heat diffusion.
     */
    buildHeatDiffusionMatrix(dt) {
        this.aDiag.fill(1);
        this.aPlusX.fill(0);
        this.aPlusY.fill(0);

        const cell = this.d.cell;
        const scale = this.diffusion * dt / (this.cellSize * this.cellSize);

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] != CELL_FLUID) continue;

                if (x < this.w - 1 && cell[idx + 1] == CELL_FLUID) {
                    this.aDiag[idx] += scale;
                    this.aDiag[idx + 1] += scale;
                    this.aPlusX[idx] = -scale;
                }
                if (y < this.h - 1 && cell[idx + this.w] == CELL_FLUID) {
                    this.aDiag[idx] += scale;
                    this.aDiag[idx + this.w] += scale;
                    this.aPlusY[idx] = -scale;
                }
            }
        }
    }

    buildPreconditioner() {
        const tau = 0.97;
        const sigma = 0.25;
        const cell = this.d.cell;

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] != CELL_FLUID) continue;

                let e = this.aDiag[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                    const px = this.aPlusX[idx - 1] * this.precon[idx - 1];
                    const py = this.aPlusY[idx - 1] * this.precon[idx - 1];
                    e = e - (px * px + tau * px * py);
                }
                if (y > 0 && cell[idx - this.w] == CELL_FLUID) {
                    const px = this.aPlusX[idx - this.w] * this.precon[idx - this.w];
                    const py = this.aPlusY[idx - this.w] * this.precon[idx - this.w];
                    e = e - (py * py + tau * px * py);
                }

                if (e < sigma * this.aDiag[idx]) {
                    e = this.aDiag[idx];
                }
                if (e < 1e-6) {
                    e = 1e-6;
                }
                this.precon[idx] = 1 / Math.sqrt(e);
            }
        }
    }

    applyPreconditioner(dst, a) {
        // applies preconditioner to vector 'a' and stores it in 'dst'
        const cell = this.d.cell;

        // first Lq = r  
        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] != CELL_FLUID) continue;

                let t = a[idx];

                if (x > 0 && cell[idx - 1] == CELL_FLUID) {
                    t -= this.aPlusX[idx - 1] * this.precon[idx - 1] * dst[idx - 1];
                }
                if (y > 0 && cell[idx - this.w] == CELL_FLUID) {
                    t -= this.aPlusY[idx - this.w] * this.precon[idx - this.w] * dst[idx - this.w];
                }

                dst[idx] = t * this.precon[idx];
            }
        }

        // solve L^T z = q
        for (let y = this.h - 1, idx = this.w * this.h - 1; y >= 0; y--) {
            for (let x = this.w - 1; x >= 0; x--, idx--) {
                if (cell[idx] != CELL_FLUID) continue;

                let t = dst[idx];

                if (x < this.w - 1 && cell[idx + 1] == CELL_FLUID) {
                    t -= this.aPlusX[idx] * this.precon[idx] * dst[idx + 1];
                }
                if (y < this.h - 1 && cell[idx + this.w] == CELL_FLUID) {
                    t -= this.aPlusY[idx] * this.precon[idx] * dst[idx + this.w];
                }

                dst[idx] = t * this.precon[idx];
            }
        }
    }

    applyA(dst, b) {
        // multipliers internal pressure matrix (A) with vecotr 'b' and stores the result in 'dst'
        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                let t = this.aDiag[idx] * b[idx];

                if (x > 0) {
                    t += this.aPlusX[idx - 1] * b[idx - 1];
                }
                if (x < this.w - 1) {
                    t += this.aPlusX[idx] * b[idx + 1];
                }
                if (y > 0) {
                    t += this.aPlusY[idx - this.w] * b[idx - this.w];
                }
                if (y < this.h - 1) {
                    t += this.aPlusY[idx] * b[idx + this.w];
                }
                dst[idx] = t;
            }
        }
    }

    projectGaussSeidel(limit, dt) {
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

    projectConjugateGradient(limit) {
        const cell = this.d.cell;

        this.pressure.fill(0);
        this.applyPreconditioner(this.z, this.rhs);
        this.s.set(this.z);

        let maxError = infinityNorm(this.rhs, cell);
        if (maxError < 1e-5) {
            console.log("Exiting solver, no pressure needed");
            return;
        }

        let sigma = dotProduct(this.z, this.rhs, cell);

        for (let iter = 0; iter < limit; iter++) {
            this.applyA(this.z, this.s);
            let alpha = sigma / dotProduct(this.z, this.s, cell);

            scaledAdd(this.pressure, this.pressure, this.s, alpha, cell);
            scaledAdd(this.rhs, this.rhs, this.z, -alpha, cell);

            maxError = infinityNorm(this.rhs, cell);
            if (maxError < 1e-5) {
                console.log(`Exiting solver after ${iter} iterations. Maximum error is ${maxError}`);
                return;
            }

            this.applyPreconditioner(this.z, this.rhs);

            const sigmaNew = dotProduct(this.z, this.rhs, cell);
            scaledAdd(this.s, this.z, this.s, sigmaNew / sigma, cell);
            sigma = sigmaNew;
        }

        console.log(`EXCEEDED MAXIMUM ITERATIONS. Maximum error is ${maxError}`);
        console.log("pres", this.pressure)
    }

    applyPressure(dt) {
        const scale = dt / this.cellSize;
        const cell = this.d.cell;

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] != CELL_FLUID) continue;

                this.u.src[this.u.id(x, y)] -= scale * this.pressure[idx] / this.uDensity[this.u.id(x, y)];
                this.u.src[this.u.id(x + 1, y)] += scale * this.pressure[idx] / this.uDensity[this.u.id(x + 1, y)];

                this.v.src[this.v.id(x, y)] -= scale * this.pressure[idx] / this.vDensity[this.v.id(x, y)];
                this.v.src[this.v.id(x, y + 1)] += scale * this.pressure[idx] / this.vDensity[this.v.id(x, y + 1)];
            }
        }
    }

    addBuoyancy(dt) {
        const alpha = (this.densitySoot - this.densityAir) / this.densityAir;

        for (let y = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++) {
                // higher density (more soot) causes downward force,
                // high temperatures cuase upward force
                const buoyancy = dt * (alpha * this.d.at(x, y) - (this.t.at(x, y) - this.tAmb) / this.tAmb);

                this.u.src[this.u.id(x, y)] += buoyancy * this.gravityX * 0.5;
                this.u.src[this.u.id(x + 1, y)] += buoyancy * this.gravityX * 0.5;

                this.v.src[this.v.id(x, y)] += buoyancy * this.gravityY * 0.5;
                this.v.src[this.v.id(x, y + 1)] += buoyancy * this.gravityY * 0.5;
            }
        }
    }

    /**
     * set all vel cells bordering solid cells to the solid velocity
     */
    setBoundaryCondition() {
        const cell = this.d.cell;
        const body = this.d.body;

        for (let y = 0, idx = 0; y < this.h; y++) {
            for (let x = 0; x < this.w; x++, idx++) {
                if (cell[idx] == CELL_SOLID) {
                    const b = this.bodies[body[idx]];

                    this.u.src[this.u.id(x, y)] = b.velocityX(x * this.cellSize, (y + 0.5) * this.cellSize);
                    this.v.src[this.v.id(x, y)] = b.velocityY((x + 0.5) * this.cellSize, y * this.cellSize);
                    this.u.src[this.u.id(x + 1, y)] = b.velocityX((x + 1) * this.cellSize, (y + 0.5) * this.cellSize);
                    this.v.src[this.v.id(x, y + 1)] = b.velocityY((x + 0.5) * this.cellSize, (y + 1) * this.cellSize);
                }
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

                const ink = this.d.at(x, y) * 100;
                const solids = (1 - this.d.volumeAt(x, y)) * 200;

                ctx.fillStyle = this.gs.getCellColour(ink - solids, 0, solids);

                /* let lkj = this.ink.phi[(this.w + 1) * y + x] * 1000;

                if (lkj <= 0){
                    ctx.fillStyle = this.gs.getCellColour(-lkj, 0, 0);
                }else{
                    ctx.fillStyle = this.gs.getCellColour(0, 0, lkj);
                } */

                /*                 const temp = this.t.at(x, y) - this.tAmb;
                                const r = temp9.81
                
                                ctx.fillStyle = this.gs.getCellColour(r, 0, 0); */

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
