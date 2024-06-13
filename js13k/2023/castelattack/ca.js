/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 p2.js authors
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


class vec2 {
	static create() { return new Float32Array([0,0]); };
	static xLen(a,b){ return a[0] * b[1] - a[1] * b[0];};
	static rotate(out,a,angle){
		var c = Math.cos(angle);
		var s = Math.sin(angle);
		var [x,y] = a;
		out[0] = c*x - s*y;
		out[1] = s*x + c*y;
	};
	static rotate90cw(out, a) { out[0] = a[1]; out[1] = -a[0];};
	static centroid(out, a, b, c){
		vec2.add(out, a, b);
		vec2.add(out, out, c);
		vec2.scale(out, out, 1/3);
	};

	static fromValues(x, y) { return new Float32Array([x, y]); };
	static copy(out, a) { out.set(a); };
	static set(out, x, y) { out.set([x, y]); };
	static add(out, a, b) { out[0] = a[0] + b[0];out[1] = a[1] + b[1];};
	static sub(out, a, b) { out[0] = a[0] - b[0]; out[1] = a[1] - b[1]; };
	static mult(out, a, b) { out[0] = a[0] * b[0];out[1] = a[1] * b[1]; };
	static scale(out, a, b) { out[0] = a[0] * b;out[1] = a[1] * b;return out;};
	static length(a) { var [x,y] = a; return Math.sqrt(x*x + y*y); };
	static sqLen(a) { var [x,y] = a;return x*x + y*y; };
	static norm(out, a) {
		var [x,y] = a;
		var len = x*x + y*y;
		if (len > 0) {
			len = 1 / Math.sqrt(len);
			out[0] = a[0] * len;
			out[1] = a[1] * len;
		}
	};
	static dot(a, b) {
		return a[0] * b[0] + a[1] * b[1];
	};
}
class AABB {
	constructor(){
		this.lowerBound = vec2.create();
		this.upperBound = vec2.create();
	}
	
	/**
	 * Set the AABB bounds from a set of points, transformed by the given position and angle.
	 * @method setFromPoints
	 * @param {Array} points An array of vec2's.
	 * @param {Array} position
	 * @param {number} angle
	 */
	setFromPoints(points, position, angle){
		var l = this.lowerBound;
		var u = this.upperBound;
		
		// Set to the first point
		vec2.rotate(l, points[0], angle);
		vec2.copy(u, l);
		
		// Compute cosines and sines just once
		var cosAngle = Math.cos(angle);
		var sinAngle = Math.sin(angle);
		for(var i = 1; i<points.length; i++){
			var [x,y] = points[i];
			var p = [
				cosAngle * x -sinAngle * y,
				sinAngle * x +cosAngle * y
			];
			for(var j=0; j<2; j++){
				if(p[j] > u[j]){
					u[j] = p[j];
				}
				if(p[j] < l[j]){
					l[j] = p[j];
				}
			}
		}
		// Add offset
		if(position){
			vec2.add(this.lowerBound, this.lowerBound, position);
			vec2.add(this.upperBound, this.upperBound, position);
		}
	};
	
	/**
	 * Copy bounds from an AABB to this AABB
	 * @method copy
	 * @param  {AABB} aabb
	 */
	copy(aabb){
		vec2.copy(this.lowerBound, aabb.lowerBound);
		vec2.copy(this.upperBound, aabb.upperBound);
	};
	
	
	/**
	 * Returns true if the given AABB overlaps this AABB.
	 * @method overlaps
	 * @param  {AABB} aabb
	 * @return {Boolean}
	 */
	overlaps(aabb){
		var l1 = this.lowerBound;
		var u1 = this.upperBound;
		var l2 = aabb.lowerBound;
		var u2 = aabb.upperBound;
		
		return ((l2[0] <= u1[0] && u1[0] <= u2[0]) || (l1[0] <= u2[0] && u2[0] <= u1[0])) &&
			((l2[1] <= u1[1] && u1[1] <= u2[1]) || (l1[1] <= u2[1] && u2[1] <= u1[1]));
	};
}


class Shape {
	constructor(options={}){
		this.body = null;
		this.position = vec2.fromValues(0,0);
		this.angle = 0;
		this.type = options.type;
		this.id = ShapeIDCounter++;
		this.material = null;
	}
}
class Box extends Shape {
	constructor(options){
		
		var width = options.width;
		var height = options.height;
		switch (options.type) {
			case "box":
			case "triangleD":
				var shape = [ //BOX
					[-0.5,-0.5],
					[0.5,-0.5],
					[0.5,0.5],
					[-0.5,0.5]
				];
				break;
			case "triangle3":
				var shape = [ //Triangle
					[-0.5,-0.5],
					[0,-0.5],
					[0.5,-0.5],
					[0,0.5]
				];
				break;
			case "triangleA":
				var shape = [ //Triangle
					[-0.5,-0.5],
					[0,-0.5],
					[0.5,-0.5],
					[0,0.]
				];
				break;
			case "triangleB":
				var shape = [ //Triangle
					[-0.5,-0.5],
					[0,-0.5],
					[0.5,-0.5],
					[0.5,0.5]
				];
				break;
			case "triangleC":
				var shape = [ //Triangle
					[-0.5,-0.5],
					[0,-0.5],
					[0.5,-0.5],
					[-0.5,0.5]
				];
				break;
				
		}
		
		options.type = 8;
		super(options);
		
		
		this.vertices = shape.map(elm => vec2.fromValues(width * elm[0], height * elm[1]));
		this.width = width;
		this.height = height;
	}
	
	cMI(mass){
		return mass * (this.height**2 + this.width**2) / 12;
	};
	computeAABB(out, position, angle){
		out.setFromPoints(this.vertices,position,angle);
	};
	
}

class Plane extends Shape {
	constructor(){
		super({type:4});
	}
	
	computeAABB(out, position, angle){
		var a = angle % (2 * Math.PI);
		
		// Set max bounds
		vec2.set(out.lowerBound, -1e7, -1e7);
		vec2.set(out.upperBound,  1e7,  1e7);
		
		if(a === 0){
			out.upperBound[1] = 0;
		} else if(a === Math.PI){
			out.lowerBound[1] = 0;
		} else if(a === Math.PI / 2){
			out.lowerBound[0] = 0;
		} else if(a === 3*Math.PI/2){
			out.upperBound[0] = 0;
		}
	}
}


class Circ extends Shape {	
	constructor(options={}){
		options.type = 1;
		super(options);		
		this.radius = options.radius;
	}
	
	cMI(mass){
		return mass * this.radius**2 / 2;
	};
	
	computeAABB(out, position, angle){
		vec2.set(out.upperBound,  this.radius,  this.radius);
		vec2.set(out.lowerBound, -this.radius, -this.radius);
		if(position){
			vec2.add(out.lowerBound, out.lowerBound, position);
			vec2.add(out.upperBound, out.upperBound, position);
		}
	};
}


class ContactMaterial {
	constructor(materialA, materialB, options={}){
		this.id = ContactMaterialIDCounter++;
		this.materialA = materialA;
		this.materialB = materialB;
		this.friction = options.friction || 0.3;
		this.restitution = options.restitution || 0;
		this.stiffness = options.stiffness || 1e6;
	}
}


class TupleDict {
	constructor() {
		this.data = {};
		this.keys = [];
	}
	getKey(id1, id2) {
		return (id1 > id2 ? (id1 << 16) | (id2 & 0xFFFF) : (id2 << 16) | (id1 & 0xFFFF))|0;
	};
	get(i, j) {
		return this.data[this.getKey(i, j)];
	};
	
	set(i, j, value) {
		var key = this.getKey(i, j);
		if(!this.data[key]){
			this.keys.push(key);
		}
		this.data[key] = value;
		return key;
	};
	
	reset() {
		var keys = this.keys;
		var l = keys.length;
		while(l--) {
			delete this.data[keys[l]];
		}
		
		keys.length = 0;
	};
	
	copy(dict) {
		this.reset();
		this.keys.push(...dict.keys);
		var l = dict.keys.length;
		while(l--){
			var key = dict.keys[l];
			this.data[key] = dict.data[key];
		}
	};
}

class OverlapKeeperRecord {
	constructor(bodyA, shapeA, bodyB, shapeB){
		this.set(...arguments);
	}
	set(bodyA, shapeA, bodyB, shapeB){
		this.shapeA = shapeA;
		this.shapeB = shapeB;
		this.bodyA = bodyA;
		this.bodyB = bodyB;
	};
}

class Equation {
	constructor(bodyA, bodyB){
		this.minForce = 0;
		this.maxForce = Number.MAX_VALUE;
		this.bodyA = bodyA;
		this.bodyB = bodyB;
		this.stiffness = 1e6;
		this.G = new Float32Array(6).fill(0);
		this.a = 0;
		this.b = 0;
		this.epsilon = 0;
		this.timeStep = 1/60;
		this.multiplier = 0;
	}
	
	update(){
		var h = this.timeStep;
		
		this.a = 4 / (h * 17);
		this.b = 16 / 17;
		this.epsilon = 4 / (h**2 * this.stiffness * 17);
	};
	
	gmult(G,vi,wi,vj,wj){
		return  G[0] * vi[0] +
			G[1] * vi[1] +
			G[2] * wi +
			G[3] * vj[0] +
			G[4] * vj[1] +
			G[5] * wj;
	};
	
	computeGW(){
		var bi = this.bodyA;
		var bj = this.bodyB;
		return this.gmult(this.G, bi.vel, bi.angVel, bj.vel, bj.angVel);
	};
	
	computeGiMf(){
		var bi = this.bodyA;
		var bj = this.bodyB;

		var iMfi = vec2.scale(vec2.create(), bi.force, bi.invMassSolve);
		var iMfj = vec2.scale(vec2.create(), bj.force, bj.invMassSolve);
		vec2.mult(iMfi, bi.massMultiplier, iMfi);
		vec2.mult(iMfj, bj.massMultiplier, iMfj);
		return this.gmult(this.G, iMfi, 0, iMfj, 0);
	};
	
	addToWlambda(deltalambda){
		var bi = this.bodyA;
		var bj = this.bodyB;
		var G = this.G;
		var Gi = vec2.fromValues(G[0], G[1]);
		var Gj = vec2.fromValues(G[3], G[4]);
		var temp = vec2.create();
		vec2.scale(temp, Gi, bi.invMassSolve*deltalambda);
		vec2.mult(temp, temp, bi.massMultiplier);
		vec2.add( bi.vlambda, bi.vlambda, temp);
		bi.wlambda += bi.invInertiaSolve * G[2] * deltalambda;
		
		vec2.scale(temp, Gj, bj.invMassSolve*deltalambda);
		vec2.mult(temp, temp, bj.massMultiplier);
		vec2.add( bj.vlambda, bj.vlambda, temp);
		bj.wlambda += bj.invInertiaSolve * G[5] * deltalambda;
	};
	
	computeInvC(eps){
		var bi = this.bodyA;
		var bj = this.bodyB;
		var invMassi = bi.invMassSolve;
		var invMassj = bj.invMassSolve;
		var invIi = bi.invInertiaSolve;
		var invIj = bj.invInertiaSolve;
		var G = this.G;
		
		return 1.0 / ((G[0] * G[0] * invMassi * bi.massMultiplier[0] +
			G[1] * G[1] * invMassi * bi.massMultiplier[1] +
			G[2] * G[2] *    invIi +
			G[3] * G[3] * invMassj * bj.massMultiplier[0] +
			G[4] * G[4] * invMassj * bj.massMultiplier[1] +
			G[5] * G[5] *    invIj) + eps);
	};
}
class ContactEquation extends Equation {
	constructor(bodyA, bodyB){
		super(bodyA, bodyB);
		
		this.conPtA = vec2.create();
		this.penetrationVec = vec2.create();
		this.conPtB = vec2.create();
		this.normalA = vec2.create();
		this.restitution = 0;
		this.firstImpact = false;
	}

	computeB(a,b,h){
		var bi = this.bodyA;
		var bj = this.bodyB;
		
		var penetrationVec = this.penetrationVec;
		var n = this.normalA;
		
		var rixn = vec2.xLen(this.conPtA,n);
		var rjxn = vec2.xLen(this.conPtB,n);
		
		this.G.set([-n[0], -n[1], -rixn, n[0], n[1], rjxn]);
		
		vec2.add(penetrationVec,bj.position,this.conPtB);
		vec2.sub(penetrationVec,penetrationVec,bi.position);
		vec2.sub(penetrationVec,penetrationVec,this.conPtA);
		
		// Compute iteration
		var GW = this.computeGW();
		var Gq = 0;
		if(this.firstImpact && this.restitution !== 0){
			GW *= (1/b)*(1+this.restitution);
		} else {
			Gq = vec2.dot(n,penetrationVec); // + 0.05;
		}
		
		return - Gq * a - GW * b - h*this.computeGiMf();
	};
}

class FrictionEquation extends Equation {
	constructor(){
		super();
		
		this.conPtA = vec2.create();
		this.conPtB = vec2.create();
		this.t = vec2.create();
		this.conEqs = [];
		this.frictionCoefficient = 0.3;
	}
	
	
	
	computeB(a,b,h){
		var t = this.t;
		
		this.G.set([-t[0], -t[1], -vec2.xLen(this.conPtA, t), t[0], t[1], vec2.xLen(this.conPtB, t)]);
		return -this.computeGW() * b - h*this.computeGiMf();
	};
}

class Pool {
	constructor(type, size) {
		this.objects = Array.from({ length: size }, () => new type);
	}
	get(type) {
		return this.objects.length ? this.objects.pop() : new type;
	};
	release(object) {
		object.bodyA = object.bodyB = null;
		this.objects.push(object);
	};
}
class OverlapKeeper {
	/**
	 * Keeps track of overlaps in the current state and the last step state.
	 * @class OverlapKeeper
	 * @constructor
	 */
	constructor() {
		this.overlappingShapesLastState = new TupleDict();
		this.overlappingShapesCurrentState = new TupleDict();
		this.recordPool = new Pool(OverlapKeeperRecord, 16);
	}
	
	/**
	 * Ticks one step forward in time. This will move the current overlap state to the "old" overlap state, and create a new one as current.
	 * @method tick
	 */
	tick() {
		var last = this.overlappingShapesLastState;
		var current = this.overlappingShapesCurrentState;
		
		// Save old objects into pool
		var l = last.keys.length;
		while(l--){
			var lastObject = last.data[last.keys[l]];
			if (lastObject){
				this.recordPool.release(lastObject);
			}
		}
		
		last.reset();
		last.copy(current);
		current.reset();
	};
	
	setOverlapping(bodyA, shapeA, bodyB, shapeB) {
		var current = this.overlappingShapesCurrentState;
		
		// Store current contact state
		if(!current.get(shapeA.id, shapeB.id)){
			var data = this.recordPool.get(OverlapKeeperRecord);
			data.set(bodyA, shapeA, bodyB, shapeB);
			current.set(shapeA.id, shapeB.id, data);
		}
	};
}



class Narrowphase {
	constructor(){
		this.conEqs = [];
		this.conEqPool = new Pool(ContactEquation, 32);
		this.fricEqs = [];
		this.frictionCoefficient = 0.3;
		this.frictionEquationPool = new Pool(FrictionEquation, 64);
		this.restitution = 0;
		this.stiffness = 1e6;
		this.collidePrev = new TupleDict();
	}
	
	reset(){
		this.collidePrev.reset();
		
		var eqs = this.conEqs;
		var l = eqs.length;
		while(l--){
			var eq = eqs[l];
			this.collidePrev.set(eq.bodyA.id, eq.bodyB.id, true);
		}
		
		for(var i=0; i<this.conEqs.length; i++){
			this.conEqPool.release(this.conEqs[i]);
		}
		for (var i=0; i<this.fricEqs.length; i++){
			this.frictionEquationPool.release(this.fricEqs[i]);
		}
		
		this.conEqs.length = this.fricEqs.length = 0;
	};
	
	createContactEquation(bodyA, bodyB, shapeA, shapeB){
		var c = this.conEqPool.get(ContactEquation);
		c.bodyA = bodyA;
		c.bodyB = bodyB;
		c.restitution = this.restitution;
		c.firstImpact = !this.collidePrev.get(bodyA.id, bodyB.id);
		c.stiffness = this.stiffness;
		return c;
	};
	
	createFrictionEquation(bodyA, bodyB, shapeA, shapeB){
		var c = this.frictionEquationPool.get(FrictionEquation);
		c.bodyA = bodyA;
		c.bodyB = bodyB;
		c.maxForce = this.slipForce;
		c.minForce = -this.slipForce;
		c.frictionCoefficient = this.frictionCoefficient;
		c.conEqs.length = 0;
		c.stiffness = 1e6;
		return c;
	};
	
	createFrictionFromContact(c){
		var eq = this.createFrictionEquation(c.bodyA, c.bodyB, c.shapeA, c.shapeB);
		
		eq.conPtA.set(c.conPtA);
		eq.conPtB.set(c.conPtB);
		
		vec2.rotate90cw(eq.t, c.normalA);
		eq.conEqs.push(c);
		return eq;
	};
	
	createFrictionFromAverage(numContacts){
		var c = this.conEqs[this.conEqs.length - 1];
		var eq = this.createFrictionEquation(c.bodyA, c.bodyB, c.shapeA, c.shapeB);
		var bodyA = c.bodyA;
		vec2.set(eq.conPtA, 0, 0);
		vec2.set(eq.conPtB, 0, 0);
		vec2.set(eq.t, 0, 0);
		for(var i=0; i!==numContacts; i++){
			c = this.conEqs[this.conEqs.length - 1 - i];
			if(c.bodyA === bodyA){
				vec2.add(eq.t, eq.t, c.normalA);
				vec2.add(eq.conPtA, eq.conPtA, c.conPtA);
				vec2.add(eq.conPtB, eq.conPtB, c.conPtB);
			} else {
				vec2.sub(eq.t, eq.t, c.normalA);
				vec2.add(eq.conPtA, eq.conPtA, c.conPtB);
				vec2.add(eq.conPtB, eq.conPtB, c.conPtA);
			}
			eq.conEqs.push(c);
		}
		var invNumContacts = 1/numContacts;
		vec2.scale(eq.conPtA, eq.conPtA, invNumContacts);
		vec2.scale(eq.conPtB, eq.conPtB, invNumContacts);
		vec2.norm(eq.t, eq.t);
		vec2.rotate90cw(eq.t, eq.t);
		return eq;
	};
	
	circCvx(circBody, circShape, circOffset, circAngle, cvxBody, cvxShape, cvxOffset,	cvxAngle){
		var radius = circShape.radius;
		var found = false;
		var minCandidateDistance = Number.MAX_VALUE;
		var verts = cvxShape.vertices;
		
		// Check all edges first
		for(var i=0; i!==verts.length+1; i++){
			var v0 = verts[i%verts.length];
			var v1 = verts[(i+1)%verts.length];
			 
			vec2.rotate(tmp1, v0, cvxAngle);
			vec2.rotate(tmp2, v1, cvxAngle);
			vec2.add(tmp1, tmp1, cvxOffset);
			vec2.add(tmp2, tmp2, cvxOffset);
			vec2.sub(tmp3, tmp2, tmp1);
			
			vec2.norm(tmp4, tmp3);
			
			// Get tangent to the edge. Points out of the Cvx
			vec2.rotate90cw(tmp5, tmp4);
			
			// Get point on circ, closest to the polygon
			vec2.scale(tmp14,tmp5,-circShape.radius);
			vec2.add(tmp14,tmp14,circOffset);
			
			if(Narrowphase.pointInCvx(tmp14,cvxShape,cvxOffset,cvxAngle)){
				vec2.sub(tmp15,tmp1,tmp14);
				var candidateDistance = Math.abs(vec2.dot(tmp15,tmp5));
				
				if(candidateDistance < minCandidateDistance){
					vec2.copy(tmp16,tmp14);
					minCandidateDistance = candidateDistance;
					vec2.scale(tmp13,tmp5,candidateDistance);
					vec2.add(tmp13,tmp13,tmp14);
					found = true;
				}
			}
		}
		
		if(found){
			var c = this.createContactEquation(circBody,cvxBody,circShape,cvxShape);
			vec2.sub(c.normalA, tmp16, circOffset);
			vec2.norm(c.normalA, c.normalA);
			
			vec2.scale(c.conPtA,  c.normalA, radius);
			vec2.add(c.conPtA, c.conPtA, circOffset);
			vec2.sub(c.conPtA, c.conPtA, circBody.position);
			
			vec2.sub(c.conPtB, tmp13, cvxOffset);
			vec2.add(c.conPtB, c.conPtB, cvxOffset);
			vec2.sub(c.conPtB, c.conPtB, cvxBody.position);
			
			this.conEqs.push(c);
			this.fricEqs.push( this.createFrictionFromContact(c) );
			return 1;
		}
		
		// Check all vertices
		for(var i=0; i<verts.length; i++){
			var localVertex = verts[i];
			vec2.rotate(tmp11, localVertex, cvxAngle);
			vec2.add(tmp11, tmp11, cvxOffset);

			vec2.sub(tmp10, tmp11, circOffset);
			if(vec2.sqLen(tmp10) < Math.pow(radius, 2)){
				var c = this.createContactEquation(circBody,cvxBody,circShape,cvxShape);

				vec2.copy(c.normalA, tmp10);
				vec2.norm(c.normalA,c.normalA);

				// Vector from circ to contact point is the normal times the circ radius
				vec2.scale(c.conPtA, c.normalA, radius);
				vec2.add(c.conPtA, c.conPtA, circOffset);
				vec2.sub(c.conPtA, c.conPtA, circBody.position);

				vec2.sub(c.conPtB, tmp11, cvxOffset);
				vec2.add(c.conPtB, c.conPtB, cvxOffset);
				vec2.sub(c.conPtB, c.conPtB, cvxBody.position);

				this.conEqs.push(c);
				this.fricEqs.push(this.createFrictionFromContact(c));
				return 1;
			}
		}
		return 0;
	};

	
	/*
	 * Check if a point is in a polygon
	 */
	static pointInCvx(worldPoint,cvxShape,cvxOffset,cvxAngle){
		var point = worldPoint;
		var verts = cvxShape.vertices;
		var lastCross = null;
		for(var i=0; i!==verts.length+1; i++){
			var v0 = verts[i%verts.length];
			var v1 = verts[(i+1)%verts.length];
			
			// Transform vertices to world
			// @todo The point should be transformed to local coordinates in the cvx, no need to transform each vertex
			vec2.rotate(tmp17, v0, cvxAngle);
			vec2.rotate(tmp18, v1, cvxAngle);
			vec2.add(tmp17, tmp17, cvxOffset);
			vec2.add(tmp18, tmp18, cvxOffset);
			
			vec2.sub(tmp19, tmp17, point);
			vec2.sub(tmp20, tmp18, point);
			var cross = vec2.xLen(tmp19,tmp20);
			
			if(lastCross===null){
				lastCross = cross;
			}
			
			// If we got a different sign of the distance vector, the point is out of the polygon
			if(cross*lastCross <= 0){
				return false;
			}
			lastCross = cross;
		}
		return true;
	}
	
	
	/**
	 * Plane/Cvx Narrowphase
	 * @method planeCvx
	 * @param  {Body} planeBody
	 * @param  {Plane} planeShape
	 * @param  {Array} planeOffset
	 * @param  {Number} planeAngle
	 * @param  {Body} cvxBody
	 * @param  {Cvx} cvxShape
	 * @param  {Array} cvxOffset
	 * @param  {Number} cvxAngle
	 */
	planeCvx(planeBody, planeShape, planeOffset, planeAngle, cvxBody, cvxShape, cvxOffset, cvxAngle){
		var numReported = 0;
		vec2.rotate(tmp2, tmp12, planeAngle);
		
		for(var i=0; i!==cvxShape.vertices.length; i++){
			var v = cvxShape.vertices[i];
			vec2.rotate(tmp1, v, cvxAngle);
			vec2.add(tmp1, tmp1, cvxOffset);
			vec2.sub(tmp3, tmp1, planeOffset);
			
			if(vec2.dot(tmp3,tmp2) <= 0){
				// Found vertex
				numReported++;
				
				var c = this.createContactEquation(planeBody,cvxBody,planeShape,cvxShape);
				vec2.sub(tmp3, tmp1, planeOffset);
				vec2.copy(c.normalA, tmp2);
				vec2.scale(tmp3, c.normalA, vec2.dot(tmp3, c.normalA));
				
				// rj is from cvx center to contact
				vec2.sub(c.conPtB, tmp1, cvxBody.position);
				
				// ri is from plane center to contact
				vec2.sub( c.conPtA, tmp1, tmp3);
				vec2.sub( c.conPtA, c.conPtA, planeBody.position);
				this.conEqs.push(c);
			}
		}
		if(numReported){
			this.fricEqs.push(this.createFrictionFromAverage(numReported));
		}
		return numReported;
	};
	
	circPlane(circBody,circShape,circOffset,ai, bj,sj,xj,aj){
		var planeBody = bj;
		var planeOffset = xj;
		var planeAngle = aj;
		
		vec2.sub(tmp1, circOffset, planeOffset);
		
		// World plane normal
		vec2.rotate(tmp2, tmp12, planeAngle);
		
		// Normal direction distance
		var d = vec2.dot(tmp2, tmp1);
		
		if(d > circShape.radius){
			return 0; // No overlap. Abort.
		}
		
		// Create contact
		var contact = this.createContactEquation(planeBody,circBody,sj,circShape);
		
		// ni is the plane world normal
		vec2.copy(contact.normalA, tmp2);
		
		// rj is the vector from circ center to the contact point
		vec2.scale(contact.conPtB, contact.normalA, -circShape.radius);
		vec2.add(contact.conPtB, contact.conPtB, circOffset);
		vec2.sub(contact.conPtB, contact.conPtB, circBody.position);
		
		// ri is the distance from plane center to contact.
		vec2.scale(tmp3, contact.normalA, d);
		vec2.sub(contact.conPtA, tmp1, tmp3 ); // Subtract normal distance vector from the distance vector
		vec2.add(contact.conPtA, contact.conPtA, planeOffset);
		vec2.sub(contact.conPtA, contact.conPtA, planeBody.position);
		
		this.conEqs.push(contact);
		this.fricEqs.push( this.createFrictionFromContact(contact) );
		return 1;
	}
	
	
	circCirc(bodyA, shapeA, offsetA, angleA, bodyB, shapeB, offsetB, angleB){
		vec2.sub(tmp1,offsetA,offsetB);
		var r = shapeA.radius + shapeB.radius;
		if(vec2.sqLen(tmp1) > Math.pow(r,2)){
			return 0;
		}
		
		var c = this.createContactEquation(bodyA,bodyB,shapeA,shapeB);
		vec2.sub(c.normalA, offsetB, offsetA);
		vec2.norm(c.normalA,c.normalA);
		
		vec2.scale( c.conPtA, c.normalA,  shapeA.radius);
		vec2.scale( c.conPtB, c.normalA, -shapeB.radius);
		
		vec2.add(c.conPtA, c.conPtA, offsetA);
		vec2.sub(c.conPtA, c.conPtA, bodyA.position);
		
		vec2.add(c.conPtB, c.conPtB, offsetB);
		vec2.sub(c.conPtB, c.conPtB, bodyB.position);
		
		this.conEqs.push(c);
		this.fricEqs.push(this.createFrictionFromContact(c));
		return 1;
	};
	
	cvxCvx(bi,si,xi,ai, bj,sj,xj,aj){
		var numContacts = 0;
		
		var found = Narrowphase.findSeparatingAxis(si,xi,ai,sj,xj,aj,tmp1);
		if(!found){
			return 0;
		}
		
		// Make sure the separating axis is directed from shape i to shape j
		vec2.sub(tmp8,xj,xi);
		if(vec2.dot(tmp1,tmp8) > 0){
			vec2.scale(tmp1,tmp1,-1);
		}
		
		// Find edges with normals closest to the separating axis
		var closestEdge1 = Narrowphase.getClosestEdge(si,ai,tmp1,true); // Flipped axis
		var closestEdge2 = Narrowphase.getClosestEdge(sj,aj,tmp1);
		
		if(closestEdge1 === -1 || closestEdge2 === -1){
			return 0;
		}
		
		// Loop over the shapes
		for(var k=0; k<2; k++){
			var closestEdgeA = closestEdge1;
			var closestEdgeB = closestEdge2;
			var shapeA = si;
			var shapeB = sj;
			var offsetA = xi;
			var offsetB = xj;
			var angleA = ai;
			var angleB = aj;
			var bodyA = bi;
			var bodyB = bj;
			
			if(k === 0){ // Swap!
				[closestEdgeA, closestEdgeB] = [closestEdgeB, closestEdgeA];
				[shapeA, shapeB] = [shapeB, shapeA];
				[offsetA, offsetB] = [offsetB, offsetA];
				[angleA, angleB] = [angleB, angleA];
				[bodyA, bodyB] = [bodyB, bodyA];
			}
			
			// Loop over 2 points in cvx B
			for(var j=closestEdgeB; j<closestEdgeB+2; j++){
				
				// Get world point
				var v = shapeB.vertices[(j+shapeB.vertices.length)%shapeB.vertices.length];
				vec2.rotate(tmp2, v, angleB);
				vec2.add(tmp2, tmp2, offsetB);
				
				var insideNumEdges = 0;
				
				// Loop over the 3 closest edges in cvx A
				for(var i=closestEdgeA-1; i<closestEdgeA+2; i++){
					
					var v0 = shapeA.vertices[(i+shapeA.vertices.length)%shapeA.vertices.length];
					var v1 = shapeA.vertices[(i+1+shapeA.vertices.length)%shapeA.vertices.length];
					
					// Construct the edge
					vec2.rotate(tmp3, v0, angleA);
					vec2.rotate(tmp4, v1, angleA);
					vec2.add(tmp3, tmp3, offsetA);
					vec2.add(tmp4, tmp4, offsetA);
					
					vec2.sub(tmp5, tmp4, tmp3);
					
					vec2.rotate90cw(tmp9, tmp5); // Normal points out of cvx 1
					vec2.norm(tmp9,tmp9);
					
					vec2.sub(tmp8, tmp2, tmp3);
					
					var d = vec2.dot(tmp9,tmp8);
					
					if((i === closestEdgeA && d <= 0) || (i !== closestEdgeA && d <= 0)){
						insideNumEdges++;
					}
				}
				
				if(insideNumEdges >= 3){
					
					// worldPoint was on the "inside" side of each of the 3 checked edges.
					// Project it to the center edge and use the projection direction as normal
					
					// Create contact
					var c = this.createContactEquation(bodyA,bodyB,shapeA,shapeB);
					numContacts++;
					
					// Get center edge from body A
					var v0 = shapeA.vertices[(closestEdgeA)   % shapeA.vertices.length];
					var v1 = shapeA.vertices[(closestEdgeA+1) % shapeA.vertices.length];
					
					// Construct the edge
					vec2.rotate(tmp3, v0, angleA);
					vec2.rotate(tmp4, v1, angleA);
					vec2.add(tmp3, tmp3, offsetA);
					vec2.add(tmp4, tmp4, offsetA);
					
					vec2.sub(tmp5, tmp4, tmp3);
					
					vec2.rotate90cw(c.normalA, tmp5); // Normal points out of cvx A
					vec2.norm(c.normalA,c.normalA);
					
					vec2.sub(tmp8, tmp2, tmp3); // From edge point to the penetrating point
					var d = vec2.dot(c.normalA,tmp8);             // Penetration
					vec2.scale(tmp7, c.normalA, d);     // Vector penetration
					
					vec2.sub(c.conPtA, tmp2, offsetA);
					vec2.sub(c.conPtA, c.conPtA, tmp7);
					vec2.add(c.conPtA, c.conPtA, offsetA);
					vec2.sub(c.conPtA, c.conPtA, bodyA.position);
					
					vec2.sub(c.conPtB, tmp2, offsetB);
					vec2.add(c.conPtB, c.conPtB, offsetB);
					vec2.sub(c.conPtB, c.conPtB, bodyB.position);
					
					this.conEqs.push(c);
				}
			}
		}
		
		if(numContacts){
			this.fricEqs.push(this.createFrictionFromAverage(numContacts));
		}
		return numContacts;
	};
	
	
	static projectCvxOntoAxis(cvxShape, cvxOffset, cvxAngle, worldAxis, result){
		var max=-Infinity;
		var min=Infinity;
		
		// Convert the axis to local coords of the body
		vec2.rotate(tmp21, worldAxis, -cvxAngle);
		
		// Get projected position of all vertices
		for(var i=0; i<cvxShape.vertices.length; i++){
			var v = cvxShape.vertices[i];
			var value = vec2.dot(v,tmp21);
			max = Math.max(max, value);
			min = Math.min(min, value);
		}
		
		if(min > max){
			[min, max] = [max, min];
		}
		
		var offset = vec2.dot(cvxOffset, worldAxis);
		vec2.set(result, min + offset, max + offset);
	};
	

	static findSeparatingAxis(c1,offset1,angle1,c2,offset2,angle2,sepAxis){
		var maxDist = null;
		var found = false;
		
		if(c1 instanceof Box && c2 instanceof Box){
			for(var j=0; j!==2; j++){
				var c = c1;
				var angle = angle1;
				if(j===1){
					c = c2;
					angle = angle2;
				}
				
				for(var i=0; i!==2; i++){
					
					// Get the world edge
					if(i === 0){
						vec2.set(tmp25, 0, 1);
					} else if(i === 1) {
						vec2.set(tmp25, 1, 0);
					}
					if(angle !== 0){
						vec2.rotate(tmp25, tmp25, angle);
					}
					
					// Project hulls onto that normal
					Narrowphase.projectCvxOntoAxis(c1,offset1,angle1,tmp25,tmp26);
					Narrowphase.projectCvxOntoAxis(c2,offset2,angle2,tmp25,tmp27);
					
					// Order by span position
					// Get separating distance
					var dist = tmp27[0] - tmp26[1];
					if(tmp26[0] > tmp27[0]){
						var dist = tmp26[0] - tmp27[1];
					}
					
					if(maxDist===null || dist > maxDist){
						vec2.copy(sepAxis, tmp25);
						maxDist = dist;
						found = dist <= 0;
					}
				}
			}
		} else {
			for(var j=0; j!==2; j++){
				var c = c1,
				angle = angle1;
				if(j===1){
					c = c2;
					angle = angle2;
				}
				
				for(var i=0; i!==c.vertices.length; i++){
					// Get the world edge
					vec2.rotate(tmp23, c.vertices[i], angle);
					vec2.rotate(tmp24, c.vertices[(i+1)%c.vertices.length], angle);
					
					sub(tmp22, tmp24, tmp23);
					
					// Get normal - just rotate 90 degrees since vertices are given in CCW
					vec2.rotate90cw(tmp25, tmp22);
					vec2.norm(tmp25,tmp25);
					
					// Project hulls onto that normal
					Narrowphase.projectCvxOntoAxis(c1,offset1,angle1,tmp25,tmp26);
					Narrowphase.projectCvxOntoAxis(c2,offset2,angle2,tmp25,tmp27);
					
					// Order by span position
					var dist = tmp27[0] - tmp26[1];
					if(tmp26[0] > tmp27[0]){
						var dist = tmp26[0] - tmp27[1];
					}
					
					if(maxDist===null || dist > maxDist){
						vec2.copy(sepAxis, tmp25);
						maxDist = dist;
						found = dist <= 0;
					}
				}
			}
		}
		return found;
	};
	
	static getClosestEdge(c,angle,axis,flip){
		// Convert the axis to local coords of the body
		vec2.rotate(tmp28, axis, -angle);
		if(flip){
			vec2.scale(tmp28,tmp28,-1);
		}
		
		var closestEdge = -1;
		var N = c.vertices.length;
		var maxDot = -1;
		for(var i=0; i!==N; i++){
			// Get the edge
			vec2.sub(tmp29, c.vertices[(i+1)%N], c.vertices[i%N]);
			
			// Get normal - just rotate 90 degrees since vertices are given in CCW
			vec2.rotate90cw(tmp30, tmp29);
			vec2.norm(tmp30,tmp30);
			
			var d = vec2.dot(tmp30,tmp28);
			if(closestEdge === -1 || d > maxDot){
				closestEdge = i % N;
				maxDot = d;
			}
		}
		return closestEdge;
	};
}

class GSSolver {
	constructor(){
		this.type = 1;
		this.lambda = new Float32Array(30);
		this.Bs =     new Float32Array(30);
		this.invCs =  new Float32Array(30);
		
	}
	
	/**
	 * Solve the system of equations
	 * @method solve
	 * @param  {Number}  h       Time step
	 * @param  {World}   world    World to solve
	 */
	solve(dt, equations){ 
		var Neq = equations.length;
		var bodies = world.bodies;
		var Nbodies = world.bodies.length;
		
		for(var i=0; i!==Nbodies; i++){
			var b = bodies[i];
			// Update solve mass
			if(bodies[i].sleepState === 2){
				b.invMassSolve = 0;
				b.invInertiaSolve = 0;
			} else {
				b.invMassSolve = b.invMass;
				b.invInertiaSolve = b.invInertia;
			}
		}
		
		// Things that does not change during iteration can be computed once
		if(this.lambda.length < Neq){
			this.lambda = new Float32Array(Neq + 30);
			this.Bs = new Float32Array(Neq + 30);
			this.invCs = new Float32Array(Neq + 30);
		}
		
		
		this.lambda.fill(0.0);
		for(var i=0; i!==equations.length; i++){
			var c = equations[i];
			c.timeStep = dt;
			c.update();
			this.Bs[i] = c.computeB(c.a,c.b,dt);
			this.invCs[i] =  c.computeInvC(c.epsilon);
		}
		
		
		if(Neq !== 0){
			for(var i=0; i!==Nbodies; i++){
				bodies[i].vlambda=[0,0];
				bodies[i].wlambda = 0;
			}
			
			// Iterate over all equations
			for(var iter=0; iter!==10; iter++){
				// Accumulate the total error for each iteration.
				var deltalambdaTot = 0.0;
				
				for(var j=0; j!==Neq; j++){
					var c = equations[j];
					
					// Compute iteration
					var lambdaj = this.lambda[j],
					GWlambda = c.gmult(c.G,c.bodyA.vlambda,c.bodyA.wlambda,c.bodyB.vlambda,c.bodyB.wlambda);
					var maxForce = c.maxForce;
					var minForce = c.minForce;

					var deltalambda = this.invCs[j] * ( this.Bs[j] - GWlambda - c.epsilon * lambdaj );

					// Clamp if we are not within the min/max interval
					if(lambdaj + deltalambda < minForce*dt){
						deltalambda = minForce*dt - lambdaj;
					} else if(lambdaj + deltalambda > maxForce*dt){
						deltalambda = maxForce*dt - lambdaj;
					}
					this.lambda[j] += deltalambda;
					c.addToWlambda(deltalambda);
					
					deltalambdaTot += Math.abs(deltalambda);
				}
				
				// If the total error is small enough - stop iterate
				if(deltalambdaTot*deltalambdaTot <= 0){
					break;
				}
			}
			
			// Add result to vel
			for(i=0; i!==Nbodies; i++){
				vec2.add( bodies[i].vel, bodies[i].vel, bodies[i].vlambda);
				bodies[i].angVel += bodies[i].wlambda;
			}
			
			var l = equations.length;
			while(l--){
				equations[l].multiplier = this.lambda[l] * 1/dt;
			}
		}
	};
}

class Body {
	constructor(options){
		this.id = ++BodyIDCounter;
		this._ID = typeof(options.ID)!="undefined" ? options.ID : null;
		if (this._ID!==null) {
			this.Y = true;
		}
		
		this.mass = options.mass || 0;
		this.invMassSolve = 0;
		this.invInertiaSolve = 0;
		
		this.massMultiplier = vec2.create();
		this.position = vec2.fromValues(options.position[0],options.position[1]);
		
		this.vel = vec2.fromValues(0,0);
		this.vlambda = vec2.fromValues(0,0);
		this.wlambda = 0;
		this.angle = options.angle || 0;
		
		this.angVel = 0;
		this.force = vec2.create();
		
		this.type = options.mass ? 1 : 2;
		this.aabb = new AABB();
		
		this.wantSleep = false;
		this.sleepState = 0;
		
		this.sleepSpeedLimit = 0.15;  //0.85 c'est pas mal aussi
		this.sleepTimeLimit = 2;
		
		this.idleTime = 0;
		this.narrowWakeUp = false;
		
		options.shape.body = this;
		this.shape = options.shape;
		
		
		// Calcul de l'inertie
		if(this.type === 2){
			this.mass = Number.MAX_VALUE;
			this.invMass = 0;
			this.inertia = Number.MAX_VALUE;
			this.invInertia = 0;
		} else {
			this.inertia = this.shape.cMI(this.mass) + this.mass*vec2.sqLen(this.shape.position);
			this.invInertia = this.inertia>0 ? 1/this.inertia : 0;
			this.invMass = 1 / this.mass;
			this.massMultiplier.set([1,1]);
		}
		
		world.bodies.push(this);
	}
	
	
	/**
	 * Called every timestep to update internal sleep timer and change sleep state if needed.
	 */
	sleepTick(dt){
		if (this.type === 2){
			return;
		}
		
		var speedSquared = vec2.sqLen(this.vel) + Math.pow(this.angVel,2);
		var speedLimitSquared = Math.pow(this.sleepSpeedLimit,2);
		
		if(speedSquared >= speedLimitSquared){
			this.idleTime = 0;
			this.sleepState = 0;
		} else {
			this.idleTime += dt;
			this.sleepState = 1;
			if(this.idleTime > this.sleepTimeLimit){
				this.sleepState = 2;
				this.angVel = 0;
				this.vel = [0,0];
				this.force = [0,0];
			}
		}
	};
	
	/**
	 * Move the body forward in time given its current vel.
	 */
	integrate(dt){ //Utilisé une seule fois, mais pas bon de le déplacer, car le this à répétition permet de gagner 15octets dans le zip
		vec2.scale(tmp31, this.force, dt * this.invMass);
		vec2.mult(tmp31, this.massMultiplier, tmp31);
		vec2.add(this.vel, tmp31, this.vel);
		
		// Regular position update
		vec2.scale(tmp32, this.vel, dt);
		vec2.add(this.position, this.position, tmp32);
		this.angle += this.angVel * dt;
		
		vec2.rotate(tmp, this.shape.position, this.angle);
		vec2.add(tmp, tmp, this.position);
		this.shape.computeAABB(shapeAABB, tmp, this.shape.angle + this.angle);
		this.aabb.copy(shapeAABB);
	};
}

class World {
	constructor(){
		this.bodies = [];
		this.solver = new GSSolver();
		this.narrowphase = new Narrowphase();
		var G = 0;
		this.gravity = vec2.fromValues(0, -G); //-0.978);
		
		this.defaultMaterial = {id:MaterialIDCounter++};
		this.defConMat = new ContactMaterial(this.defaultMaterial,this.defaultMaterial);
		
		this.conMats = [];
		this.stepping = false;
		this.removeBodies = [];
		this.overlapKeeper = new OverlapKeeper();
	}
	
	
	/**
	 * Step the physics world forward in time.
	 */
	step(dt){
		this.stepping = true;
		
		var bodies = this.bodies;
		
		var Nbodies = this.bodies.length;
		var np = this.narrowphase;
		
		this.overlapKeeper.tick();
		
		// Add gravity to bodies
		for(var i=0; i!==Nbodies; i++){
			var bi = bodies[i];
			if(bi.type !== 1 || bi.sleepState === 2){
				continue;
			} 
			vec2.scale(tmp33,this.gravity,bi.mass); // F=m*g
			vec2.add(bi.force,bi.force,tmp33);
		}
		
		
		/**
		 * Get the colliding pairs
		 */
		var result = [];
		
		// Sort the lists
		for(var i=1; i!==Nbodies; i++) {
			var bi = bodies[i];
			for(var j=i - 1;j>=0;j--) {
				if(bodies[j].aabb.lowerBound[0] <= bi.aabb.lowerBound[0]){
					break;
				}
				bodies[j+1] = bodies[j];
			}
			bodies[j+1] = bi;
		}
		
		// Look through the X list
		for(var i=0; i!==Nbodies; i++){
			var bi = bodies[i];
			
			for(var j=i+1; j<Nbodies; j++){
				var bj = bodies[j];
				
				// Bounds overlap?
				if (bj.aabb.lowerBound[0] > bi.aabb.upperBound[0]) {
					break;
				}
				
				// Cannot collide static bodies
				var cannotCollide = ((bi.type === 2 && bj.type === 2) || (bi.sleepState === 2 && bj.sleepState === 2) || ((bi.sleepState === 2 && bj.type === 2) || (bj.sleepState === 2 && bi.type === 2)));
				
				if(!cannotCollide && bi.aabb.overlaps(bj.aabb)){
					result.push(bi,bj);
				}
			}
		}
	
		// Narrowphase
		np.reset();
		for(var i=0; i!==result.length; i+=2){
			var bodyA = result[i];
			var bodyB = result[i+1];
			var mat = false;
			if(bodyA.shape.material && bodyB.shape.material){
				for(var j=0; j!==this.conMats.length; j++){
					var cm = this.conMats[j];
					if(
						(cm.materialA.id === bodyA.shape.material.id) && (cm.materialB.id === bodyB.shape.material.id) 
						||
						(cm.materialA.id === bodyB.shape.material.id) && (cm.materialB.id === bodyA.shape.material.id)){
						mat = cm;
						break;
					}
				}
			}
			this.runNarrowphase(
				np,
				bodyA,
				bodyA.shape,
				bodyA.shape.position,
				bodyA.shape.angle,
				bodyB,
				bodyB.shape,
				bodyB.shape.position,
				bodyB.shape.angle,
				mat ? mat : this.defConMat,
				-this.gravity[1]
			);
		}
		
		// Wake up bodies
		for(var i=0; i!==Nbodies; i++){
			var bi = bodies[i];
			if(bi.narrowWakeUp){
				bi.idleTime = 0;
				bi.sleepState = 0;
				bi.narrowWakeUp = false;
			}
		}
		
		if(np.conEqs.length || np.fricEqs.length){
			this.solver.solve(dt, np.conEqs.concat(np.fricEqs));
		}
		
		// Step forward
		for(var i=0; i!==Nbodies; i++){
			var bi = bodies[i];
			bi.integrate(dt); 
			bi.force = [0,0];
		}
		
		// Emit impact event
		for(var i=0; i!==np.conEqs.length; i++){
			var eq = np.conEqs[i];
			if(eq.firstImpact){
				impact({
					bodyA : eq.bodyA,
					bodyB : eq.bodyB
				});
			}
		}
		
		// Sleeping update
		for(i=0; i!==Nbodies; i++){
			var bi = bodies[i];
			bi.sleepTick(dt);
		}
		
		this.stepping = false;
		
		// Remove bodies that are scheduled for removal
		for(var i=0; i!==this.removeBodies.length; i++){
			this.removeBody(this.removeBodies[i]);
		}
	};
	
	/**
	 * Runs narrowphase for the shape pair i and j.
	 * @method runNarrowphase
	 * @param  {Narrowphase} np
	 * @param  {Body} bi
	 * @param  {Shape} si
	 * @param  {Array} xi
	 * @param  {Number} ai
	 * @param  {Body} bj
	 * @param  {Shape} sj
	 * @param  {Array} xj
	 * @param  {Number} aj
	 * @param  {Number} mu
	 */

	runNarrowphase(np,bi,si,xi,ai,bj,sj,xj,aj,cm,glen){
		// Get world position and angle of each shape
		vec2.rotate(tmp34, xi, bi.angle);
		vec2.rotate(tmp35, xj, bj.angle);
		vec2.add(tmp34, tmp34, bi.position);
		vec2.add(tmp35, tmp35, bj.position);
		var aiw = ai + bi.angle;
		var ajw = aj + bj.angle;
		np.frictionCoefficient = cm.friction;
		if(bi.type === 2){
			var reducedMass = bj.mass;
		} else if(bj.type === 2){
			var reducedMass = bi.mass;
		} else {
			var reducedMass = (bi.mass*bj.mass)/(bi.mass+bj.mass);
		}
		np.slipForce = cm.friction*glen*reducedMass;
		np.restitution = cm.restitution;
		np.frictionStiffness = 1e6;
		np.stiffness = cm.stiffness;
			
		switch (si.type | sj.type) {
			case 4|8:
			case 8|32:
				var resolver = np['planeCvx'];
				break;
			case 1|4:
				var resolver = np['circPlane'];
				break;
			case 8:
			case 8|32:	
			case 32:
				var resolver = np['cvxCvx'];
				break;
			case 1|8:
			case 1|32:
				var resolver = np['circCvx'];
				break;
			case 1|1:
				var resolver = np['circCirc'];
				break;
		}
			
		var numContacts = 0;
		if (resolver) {
//			var numFrictionBefore = np.fricEqs.length;
			if (si.type < sj.type) {
				numContacts = resolver.call(np, bi,si,tmp34,aiw, bj,sj,tmp35,ajw);
			} else {
				numContacts = resolver.call(np, bj,sj,tmp35,ajw, bi,si,tmp34,aiw);
			}
			
			if(numContacts){
				if( bi.type === 1 &&
					bi.sleepState  === 2 &&
					bj.sleepState  === 0 &&
					bj.type !== 2
					){
					var speedSquaredB = vec2.sqLen(bj.vel) + Math.pow(bj.angVel,2);
					var speedLimitSquaredB = Math.pow(bj.sleepSpeedLimit,2);
					if(speedSquaredB >= speedLimitSquaredB*2){
						bi.narrowWakeUp = true;
					}
				}
				
				if( bj.type === 1 &&
					bj.sleepState  === 2 &&
					bi.sleepState  === 0 &&
					bi.type !== 2
					){
					var speedSquaredA = vec2.sqLen(bi.vel) + Math.pow(bi.angVel,2);
					var speedLimitSquaredA = Math.pow(bi.sleepSpeedLimit,2);
					if(speedSquaredA >= speedLimitSquaredA*2){
						bj.narrowWakeUp = true;
					}
				}
				
				this.overlapKeeper.setOverlapping(bi, si, bj, sj);
			}
		}
	};
	
	
	
	/**
	 * Remove a body from the simulation. If this method is called during step(), the body removal is scheduled to after the step.
	 */
	removeBody(body){
		if (this.stepping){
			this.removeBodies.push(body);
			this.removed = true;
		} else {
			var idx = this.bodies.indexOf(body);
			if(idx!==-1){
				for (var i=idx, len=this.bodies.length-1; i < len; i++){
					this.bodies[i] = this.bodies[i + 1];
				}
				this.bodies.length = len;
			}
		}
	};
}



for (var i = 1; i <= 35; i++) {
    window['tmp'+i] = vec2.create();
}
var tmp12 = vec2.fromValues(0,1);
var ContactMaterialIDCounter = 0;
var BodyIDCounter = 0;
var ShapeIDCounter = 0;
var shapeAABB = new AABB();
var tmp = vec2.create();
var MaterialIDCounter = 0;
var ctx;
var w;
var h;
var world;
var boxBody;
var planeBody;
var camera;
var steps;
	
var balls = [];
var keys = [];
keys[38] = false;
keys[40] = false;

var areas = [{
	ammos:["BOULDER", "BOULDER", "BOULDER", "BOULDER", "BOULDER"],
	level:`
                       X
                      BDDC
                      W11E
                X     WMME
               BWE    W11E
        X     BDWExxxxW11E
       BWEC   W11E    W11E
       WW1ExxxWMME    W11EC
       WWME111W11E   BW11EE
       WW1E111W11EDDDDW11EE
       WW1ExxxW11E1111W11EE
`.split('\n')
},{
	ammos:["IMPACT BOMB", "IMPACT BOMB", "IMPACT BOMB", "IMPACT BOMB", "IMPACT BOMB"],
	level:`
          X
          WE          
          WE          
         3WE3          
         W1EE         
         WMEE         
         W1EE         
         WMEE         
         W1EE      X
         WMEE   BDDDDDDC
         W1EE  BDDDDDDDDC
        BW1EE  W11111111E
    BDDDDW1EE  W1M1451M1E
    W1111W1EE  W11167111E
    W1111W1EEX W11189111E
`.split('\n')
},{ //5
	ammos:["IMPACT BOMB",  "BOULDER",  "IMPACT BOMB", "TIME BOMB"],
	level:`
                                 X
                              W1 1 1E
                              W11111E
                              W1M1M1E
      BC      A      A        W11111E
     BDDC    BDC    BDC    BC W11111E
     W11E    W1E   BDDDC   WE W11111E
     WMME X  WME   W111E  BWECW11111E
BDDDDW11E 1 1W1Ec( WM1ME  WWEEW11111E
W111EW11E1111111111W111EX WWEEW11111E
W111EW11E1111451111W111E11WWEEW11111E
W111EW11E1111671111W111E11WWEEW11111E
W111EW11E1111891111W111E11WWEEW11111E
`.split('\n')
},{ //4
	ammos:["IMPACT BOMB", "TIME BOMB", "BOULDER",  "IMPACT BOMB", "IMPACT BOMB"],
	level:`
              3                  A
             BDC                BDC
             W1E                W1E
             WME                WME
             W1E                W1E
             W1E                W1E
            BDDDC              BDDDC
            W111E              W111E
            WM1ME              WM1ME
            W111E              W111E
           BDDDDDC            BDDDDDC
           W11111E            W11111E
           W11111E   BDDDDC   W11111E  
           W11111E X W1111EX  W11111E
   x       W11111ExxxW1451ExxxW11111E
  xxX     xW11111ExxxW1671ExxxW11111Ex
 xxxW1E xxxW11111ExxxW1891ExxxW11111Exxx
`.split('\n')
},{ //7
	ammos:["TIME BOMB",  "IMPACT BOMB",  "IMPACT BOMB",  "IMPACT BOMB",  "IMPACT BOMB"],
	 camera:{x:-160, y:190},
	 zoom:20,
	 level:`

               BC
              BDDC
             BDDDDC
             W1111E
             WM11ME
             W1111E
         BC  W1xx1E  BC
        BDDC W1  1E BDDC
        W11E W1X 1E W11E
W 1 1 E WxxE W1111E WxxE W 1 1 E
W11111E W  E W1111E W  E W11111E
WMMMMME WX E W1111E WX E WMMMMME
W11111E W11E W1111E W11E W11111E
W11111E W11E W1111E W11E W11111E
W11111ExxxxxxxxxxxxxxxxxxW11111E
W11111E111111111111111111W11111E
W11111E11M11M114511M11M11W11111E
W11111E111111116711111111W11111E 
W11111E111111118911111111W11111E     
`.split('\n')
},{ //3
	ammos:["IMPACT BOMB", "IMPACT BOMB", "TIME BOMB", "TIME BOMB", "TIME BOMB"],
	level:`
            3
           BDC
          BDDDC
          W111E
          W1M1E
 BDDDDDDDDW111E
 W11111111W1M1Exxxxxxxxx
 W1M1111M1W111Exxxxx   x
 W11111111xxxxxx     X x
 xxxxxxxxxxxxxxx     WE 
                     WE 
                     WE 
                   xxxxx
                       x     A
                       x    BDC
                xW1Ex  x    W1E
                xW1Ex  xDDDDWME
 xx   BDDDC     xWMEx  x1111W1E  
 xx   W111Exxx  xW1Ex  x1M11W1E
 xx   W1M1Exxx  xW1Ex  x1111DDDC
 xx   W111ExxxX xW1ExX x1111111E
`.split('\n')
},{ //6
	ammos:["IMPACT BOMB",  "IMPACT BOMB",  "IMPACT BOMB", "TIME BOMB"],
	level:`
                              3
                             BDC
                             W1E
                     A       WME
                    BDC      W1E  BDC 
                   BDDDC     W1E  W1E 
xxx                W1xxxxxx  xxx  xxx
xxx                WMxxx
xxx BDDDC          W1xxx
xxxxxx11E          WMxxx
xxxx  1ME        BDW1xxx                 
xxxxX 11E        W1xxxxx
xxxxxxxxxx       WME xxx 1 1 E
xxxx             W1E xxx11111E
xxxxBDC          W1E xxx1MMMME
xxxxW1EBC       BDDDCxxx11111E
xxxxWMEWE       W111Exxx11111E
xxxxW1EWE       W111Exxx11451E
xxxxW1EWEDDDC   W111Exxx11671E
xxxxW1EWEW11E X W111Exxx11891E    X  
`.split('\n')
},];



function Particle(x,y,z,velX, velY, velZ) {
	this.x = x;
	this.y = y;
	this.z = z;
	this.velX = velX;
	this.velY = velY;
	this.velZ = velZ;
	this.type = "particle";
	this.frame = 0;
	this.update = function() {
		this.frame++;
		this.velY -=0.001;
		this.x+=this.velX;
		this.y+=this.velY;
		this.z+=this.velZ;
		
		ctx.beginPath();
		var x = this.x;
		var y = Math.max(-0.5, this.y);
		
		var radius = 0.05;
		ctx.strokeStyle = 'rgba(255, 139, 50, '+(1-(1/60*this.frame))+')'; 
		ctx.arc(x*zoom,y*-zoom,radius * zoom,0,2*Math.PI);
		ctx.stroke();
		
		if (this.frame==60) {
			return false;
		}
		return true;
	};
};
function addExplosion(x,y,z,nb,power) {
	for (var i=0; i<nb; i++) {
		var velX = (Math.random() * power) - (power/2);
		var velY = (Math.random() * power) - (power/2);
		var velZ = (Math.random() * power) - (power/2);
		objects.push(new Particle(x, y, z, velX, velY, velZ));
	}
};
function addObject(obj) {
	obj.ID = objects.length;
	obj.addImpact = function(i, fromExplode) {
		this.life-=i;
		if (this.life<0) {
			objects[this.ID].destroy();
		} else if (this.type=="box") {
			zzfx(1,.05,148,.01,.01,.16,4,.49,2,0,0,0,0,.4,0,.4,0,.5,0,.06); // Hit 98
		}
		
		if (this.type=="circ" && this.weapon=="IMPACT BOMB" && !fromExplode) {
			this.explode();
		}
	};
	obj.destroy = function() {
		if (this.type=="box" && this.height==1) {
			addExplosion(this.body.position[0], this.body.position[1], 0, 10, 0.02);
			zzfx(2.44,.05,615,0,.18,.5,2,3.62,.7,0,0,0,0,1.8,-44,.5,.03,.32,.16,.22); // Explosion 87
		}
		if (!this.died && this.type=="box" && this.height==1.5) { //Chests
			addExplosion(this.body.position[0], this.body.position[1], 0, 10, 0.12);
			zzfx(1.36,0,261.6256,.02,.17,.58,0,.6,0,.2,250,.09,.05,0,0,0,.06,.5,.4,0); // Pickup 159
			nbChests--;
			if (!nbChests) {
				console.log("Je pense qu'il ne reste plus de chest");
				end();
			}
		}
		this.died = true;
		world.removeBody(this.body);
		if (this.type=="box" && this.height==1 && Math.random()>0.5) {
			var x = this.body.position[0];
			var y = this.body.position[1];
			var subType = this.subType;
			var part = this.part;
			
			addObject({
				type:"box",
				position:[x, y],
				width:0.5,
				height:0.5,
				subType,part

			});
			addObject({
				type:"box",
				position:[x+0.5, y],
				width:0.5,
				height:0.5,
				subType,part
			});

			addObject({
				type:"box",
				position:[x, y-0.5],
				width:1,
				height:0.5,
				subType,part

			});
		}
	};
	obj.explode = function() {
		var x = this.body.position[0];
		var y = this.body.position[1];
		addExplosion(x, y, 0, 100, 0.2);

		if (!explode.frame) {
			explode.position.x = camera.x;
			explode.position.y = camera.y;
		}
		explode.frame = 1;
		for (var j=0; j<world.bodies.length; j++) {
			var b = world.bodies[j];
			if (b.mass) {
				var x2 = b.position[0];
				var y2 = b.position[1];

				var angle = Math.atan2(x2-x, y2-y);
				f = 50-Math.min(Math.hypot(x-x2, y-y2)*10,50);

				b.force = [Math.sin(angle)*f, Math.cos(angle)*f];

				if (b.Y && b._ID!=this.ID) {
					objects[b._ID].addImpact(f, true);
				}
			}
		}
		world.removeBody(this.body);
		objects[this.ID].destroy();
		balls.splice(balls.indexOf(this),1);
		for(var i=0; i!==world.bodies.length; i++){
			world.bodies[i].idleTime = 0;
		}
	}
	objects.push(obj);
	obj.life = obj.life || 50;
	switch (obj.type) {
		case "triangle3": 
		case "triangleA":
		case "triangleB":
		case "triangleC":
		case "triangleD":
			obj.width = 1;
			obj.height = 1;
			obj.mass = 0.1;
			obj.shape = new Box({type:obj.type, width: obj.width, height: obj.height});
			obj.shape.material = materialCrate;
			break;
		case "box":
			obj.width = obj.width || 1;
			obj.height = obj.height || 1;
			obj.mass = obj.mass===0 ? 0 : 0.1;
			obj.shape = new Box({type:"box", width: obj.width, height: obj.height});
			obj.shape.material = materialCrate;
			break;
		case "circ":
			obj.life = 1000;
			obj.shape = new Circ({radius: obj.radius});
			obj.shape.material = materialBall;
			break;
		default:
			break;
	}
	obj.body = new Body(obj);
	return obj;
};
materialBall = {id:MaterialIDCounter++};
materialCrate= {id:MaterialIDCounter++};

function drawAll() {
	for (var i=0; i<objects.length; i++) {
		if (objects[i].died) { continue; }
		switch (objects[i].type) {
			case "box":
			case "triangle3": 
			case "triangleA":
			case "triangleB":
			case "triangleC":
			case "triangleD":
				drawBox(objects[i]);
				break;
			case "circ":
				drawCirc(objects[i]);
				break;
			case "particle":
				var state = objects[i].update();
				if (!state) {
					objects[i].died = true;
				}
				break;
		}
	}
};

function drawBox(obj){
	if (obj.shape.vertices) {
		var x = obj.body.position[0];
		var y = obj.body.position[1];
		ctx.save();
		ctx.translate(x*zoom, y*-zoom);        // Translate to the center of the box
		ctx.rotate(-obj.body.angle);  // Rotate to the box body frame
		vertices = obj.shape.vertices;

		ctx.beginPath();
		var firstVertex = vertices[0];
		ctx.moveTo(firstVertex[0]*zoom, -firstVertex[1]*zoom);

		for (var i = 1; i < vertices.length; i++) {
			ctx.lineTo(vertices[i][0]*zoom, -vertices[i][1]*zoom);
		}

		if (obj.subType=="X") { //Coffre
			ctx.translate(zoom/1, zoom/1.35);
			ctx.scale(zoom/40, zoom/40);
		} else {
			ctx.translate(-0.51*zoom, -0.51*zoom);
			ctx.scale(zoom/40, zoom/(obj.height*40));
		}
		if (obj.type.substring(0,8)=="triangle") {
			ctx.fillStyle = textures.T;
		} else if (obj.subType=="windows") {
			ctx.fillStyle = textures.w[obj.part];
		} else {
			ctx.fillStyle = textures.b[50-Math.floor(obj.life)];
			if (obj.subType) {
				ctx.fillStyle = textures[obj.subType];
			}
		}
		ctx.fill();
		ctx.restore();
	}
	
	
	if (false) { //#DEBUG
		ctx.save(); //#DEBUG
		ctx.translate(x*zoom, y*-zoom); //#DEBUG
		ctx.rotate(-obj.body.angle); //#DEBUG
//		ctx.strokeStyle = obj.sleep ? "blue" : "red";//#DEBUG
		var texte = obj.life;//#DEBUG
		ctx.font = "12px Arial";//#DEBUG
		ctx.fillStyle = "red"; //#DEBUG
		ctx.fillText(texte.toFixed(0), -obj.shape.width/4*zoom, 0);//#DEBUG
//		ctx.stroke();//#DEBUG
		ctx.restore();//#DEBUG
	} //#DEBUG
};
function drawCirc(obj){
	ctx.fillStyle = textures.R;
	var radius = obj.shape.radius*2;
	fillForm(obj.body.position[0], obj.body.position[1], radius, radius, radius);
};


/**
 * Creates a pseudo-random value generator. The seed must be an integer.
 *
 * Uses an optimized version of the Park-Miller PRNG.
 * http://www.firstpr.com.au/dsp/rand31/
 */
function random(s) {
	this.s = s % 2147483647;
	this.n = function () {
		this.s = this.s * 16807 % 2147483647;
		return (this.s - 1) / 2147483646;
	};
}


//
// Generate a mountain and draw it in a canva
//
var mountain = document.createElement("canvas");
mountain.width = 1008;
mountain.height = 480;
document.body.appendChild(mountain);
function generateMountain(s) {
	v = 480;
	t = [(mr.n() * v)|0, (mr.n() * v)|0];
	while(t.length < 1008) {
		v *= 0.52;
		t = function(o, v) {
			n = [o[0]];
			for(i = 1; i < o.length; i++) {
				n.push(((o[i - 1] + o[i]) / 2 + ((mr.n() - 0.5) * v))|0);
				n.push(o[i]);
			}
			return n;
		}(t, v);
	}
	var ctx = mountain.getContext("2d");
	var grd=ctx.createLinearGradient(0,480,0,0);
	grd.addColorStop(0,"#000");
	grd.addColorStop(0.95,"rgb(100,"+Math.floor(100+mr.n()*155)+",100)");
	ctx.strokeStyle = grd;

	for(var i = 0; i < 1008; i++) {
		ctx.beginPath();
		ctx.moveTo(i + 0.5, 480 - t[i]);
		ctx.lineTo(i + 0.5, 480);
		ctx.stroke();
	}
}
var world;	
var clouds = [];
var cloud = new Image();
cloud.src = "images/c.png";
cloud.onload = function() {
	for (var i=0; i<10; i++) {
		clouds.push({x:Math.random()*1224, y:Math.random()*250, size:0.5+Math.random()/2, speed:Math.random()/4});
	}
}
var images = [];
// Init p2.js
// Configurez la gestion des collisions
	function impact(event) {
		var shapeA = event.bodyA.shape;
		var shapeB = event.bodyB.shape;

		// Calculez la vitesse relative des objets en collision
		var relativeVel = [shapeB.body.vel[0] - shapeA.body.vel[0], shapeB.body.vel[1] - shapeA.body.vel[1]];

		// Calculez la norme de la vitesse relative pour obtenir une estimation de la force d'impact
		var impactForce = Math.sqrt(relativeVel[0] * relativeVel[0] + relativeVel[1] * relativeVel[1]);

		if (impactForce<1) {
			return;
		}
		if (shapeA.body.Y) {
			objects[shapeA.body._ID].addImpact(impactForce);
		}
		if (shapeB.body.Y) {
			objects[shapeB.body._ID].addImpact(impactForce);
		}
	};
	function initLevel(l) {
		nbChests = 3;
		level = areas[l].level;
		ammos = areas[l].ammos.slice();
		world = new World();
		objects = [];

		// Level               
		for (var x=0; x<level[level.length-2].length; x++) {
			for (var y=0; y<level.length; y++) {
				var l = level[y][x];
				switch (l) {
					case "X": //Trésor
						addObject({
							type:"box",
							subType:l,
							position:[x+0.5, level.length-y-1.75],
							width:2,
							height:1.5,
							life:10
						});
						break;
					case "1":
					case "M":
					case "E":
					case "W":
					case "x":
						addObject({
							type:"box",
							subType:l,
							position:[x, level.length-y-2],
							mass:l=="x" ? 0 : 0.1,
							life:l=="x" ? 50000 : 50
						});
						break;
					case "c":
						addObject({
							type:"box",
							subType:l,
							position:[x, -0.25+(level.length-y-2)],
							width:1,
							height:0.5,
							angle:(level.length-y-1)%2==0 ? 0 : Math.PI
							
						});
						break;
					case "(":
						addObject({
							type:"box",
							subType:l,
							position:[x+0.5, -0.25+(level.length-y-2)],
							width:1,
							height:0.5,
							angle:(level.length-y-1)%2==0 ? 0 : Math.PI
						});
						break;
					case "3":
					case "A":
					case "B":
					case "C":
					case "D":
						addObject({
							type:"triangle"+(l),
							position:[x, level.length-y-2]
						});
						break;
					case " ":
					case undefined:
						break;
					default:
						addObject({
							type:"box",
							subType:"windows",
							position:[x, level.length-y-2],
							width:1,
							height:1,
							angle:l==5 ? -Math.PI/2 : 0,
							part:l-4							
						});
						break;
				}
			}			
		}

		planeShape = new Plane();
		planeBody = new Body({position:[0,-0.5],angle:0, shape:planeShape});
		planeShape.material = {id:MaterialIDCounter++};
		world.conMats.push(new ContactMaterial(planeShape.material, materialBall, {
			restitution : 0.8,
			stiffness : 20000000,    // This makes the contact soft!
			friction:100
		}));

		world.conMats.push(new ContactMaterial(materialCrate, materialCrate, {
			restitution : 0,
			stiffness : 10000000,    // This makes the contact soft!
			friction:1
		}));

		world.conMats.push(new ContactMaterial(planeShape.material, materialCrate, {
			restitution : 0,
			stiffness : 100000,    // This makes the contact soft!
			friction:10000
		}));
	}


function init(){
	// Init canvas
	//canvas = document.getElementById("myCanvas");
	w = canvas.width;
	h = canvas.height;

	ctx = canvas.getContext("2d");

	// Création du motif de brique procédural en utilisant ImageData
	function createTextures(crackSize, type, part, width = 40, height = 40) {
		var patternC = document.createElement('canvas');
		images[type+part] = patternC;
		patternC.width = width;
		patternC.height = height;
		var patternCtx = patternC.getContext('2d');
		var iData = patternCtx.createImageData(width, height);
		var data = iData.data;
		
		var jointColor = [100, 80, 70];
		var metalColor = [50, 40, 35];
		mr = new random(1001);
		
		for (var y = 0; y < height; y++) {
			for (var x = 0; x < width; x++) {
				var pixelIndex = (y * width + x) * 4;
				var r = Math.random() * 20 - 10; // Variation de couleur entre -10 et 30
				var color = [190 + r, 150 + r, 100 + r];
				
				if ((y==height/2) || (y>height/2 && x==width/2) || (y<height/2 && x==0) || y === 0) {
					color = jointColor;
				}

				if (type=="E") { //Effet d'ombre à droite
					if (x>25) {
						percent = 1-(1/30*(x-25));
						color = [color[0]*percent, color[1]*percent, color[2]*percent];
					}
				}
				if (type=="W") { //Effet d'ombre à gauche
					if (x<15) {
						percent = 0.5+(1/30*x);
						color = [color[0]*percent, color[1]*percent, color[2]*percent];
					}
				}
				if (type=="M") { //Meurtriere
					if (x>14 && x<26) {
						percent = 0;
						color = [color[0]*percent, color[1]*percent, color[2]*percent];
					}
				}

				if (type=="tuile") {
					var isBlackTileX = (y%20<10 && (x==20 || x==0)) || (y%20>10 && (x==10 || x==30));
					var isBlackTileY = Math.floor(y) % 10 === 0;
					const colorValue = isBlackTileX === isBlackTileY ? 10+Math.floor(Math.random()*30) : 80;
					color = [colorValue,colorValue,colorValue];
				}

				/**/
				if (type=="windows") {
					// Calculate wood-like color based on noise function
					var r = Math.floor(Math.random() * 70) + 145;
					var g = Math.floor(Math.random() * 30) + 60;
					var b = Math.floor(Math.random() * 30) + 10;      

					switch (part) {
						case 1:
						case 2:
							var dist = Math.hypot(80-x, 80-y);
							if (dist<90) {
								color = [r*0.9, g*0.9, b*0.9]; //La multiplication par 0.8 est là pour simuler l'ombre derrière la grille.
								if ((x+3)%10<3 || (y+3)%10<3) {
									color = metalColor;
								}
								if (dist>87.75) {
									color = metalColor;
								} 
							} else if (x==width/2) {
								color = jointColor;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
							color = [r+x, g+x, b+x];
							if (part==4 || part==6) {
								x = 40-x;
							}

							if (part==3 || part==5 ||part==4 || part==6) {
								if (x<5) { //Dégradé gauche
									percent = 1/3*x;
									color = [color[0]*percent, color[1]*percent, color[2]*percent];
								}
							}
							if (x%10<2) { //Lames
								percent = 0.8;
								color = [color[0]*percent, color[1]*percent, color[2]*percent];
							}	

							if (part<5) {
								if (y<10) { //Pointe porte
									if ((x+5)%20>y/2 && (x+5)%20<y/2 + 10 - (y)) {
										color = metalColor;
									}
								}
							}

							if (x>37 || (y>20 && y<30)) { //Barre renfort porte //Barre verticale centrale
								color = metalColor;
							}

							if (part==5 || part==6) { //Poignées portes
								var dist = Math.hypot(30-x, 5-y);
								if (dist<5 && dist>2.5) {
									color = metalColor;
								}
							}

							if (part==4 || part==6) {
								x = 40-x;
							}
							break;
						case 7:
							color = [r*0.9, g*0.9, b*0.9];
							if (x==2 || x==237 || y==2 || y==17 || ((y>8 && y<12) && x%10<3) || (x>237 || y>17)) {
								color = [r*0.6, g*0.6, b*0.6];
							}
							if (x<2 || y<2) {
								color = [r*1.2, g*1.2, b*1.2];
							}
							break;
						case 8:
							color = [r*1.2,g*2.4,b*2.5];
							break;
						case 9: //COFFRE
							//80x60
							color = [r,g,b];
							var metal = false;
							if (x<8 || x>80-8 || y<4 || (y>18 && y<26) || y>60-8) { //Contour
								metal = true;
							}

							if ((y<18 && y%8<2) || ((y+26)%10 < 2)) { //Planche
								color = [r*0.8, g*0.8, b*0.8];
							}
							
							var dist = Math.hypot(20-x, 39-y); //Rond 1
							if (dist<3) {
								color = metalColor;
							}
							var dist = Math.hypot(60 - x, 39 - y); //Rond 2
							if (dist < 3 || (y>25 && y<35 && x>36 && x<44)) { //Serrure
								metal = true;
							}
							
							if (metal) {
								color = metalColor;
								if ((x%2==0 && y%2==0) || (x%2==1 && y%2==1)) {
									color = [color[0]*0.5, color[1]*0.5, color[2]*0.5];
								}
							}
							break;
						case 10: //ROCK
							color = [color[0]*0.4, color[1]*0.5, color[2]*0.6];
					}
				}				
				data[pixelIndex] = color[0];
				data[pixelIndex + 1] = color[1];
				data[pixelIndex + 2] = color[2];
				data[pixelIndex + 3] = 255;
			}
		}
		
		var cracks = [];
		for (var f=0; f<20; f++) {
			cracks.push([
				Math.floor(mr.n()*40),
				Math.floor(mr.n()*40)
			]);
		}
		cracks.forEach(function(crack, n) {
			mr = new random(n);
			var x = crack[0];
			var y = crack[1];
			
			for (var f=0; f<crackSize; f++) { //Longueur
				x+=mr.n()<0.33 ? -1 : (mr.n()<0.66 ? 0 : 1);
				y+=mr.n()<0.33 ? -1 : (mr.n()<0.66 ? 0 : 1);
				
				pixelIndex = (y * width + x) * 4;
				color = [152,120,80];

				data[pixelIndex] = color[0];
				data[pixelIndex + 1] = color[1];
				data[pixelIndex + 2] = color[2];
				data[pixelIndex + 3] = 255;
			}
		});
		
		patternCtx.putImageData(iData, 0, 0);
		document.body.appendChild(patternC);//#DEBUG
		return ctx.createPattern(patternC, 'repeat');
	}

	textures = {
		X:createTextures(0, "windows", 9, 80,60), //Trésor
		T:createTextures(0, "tuile", 1), //Tuile
		W:createTextures(0, "W", 1), //Ombre à gauche
		E:createTextures(0, "E", 1), //Ombre à droite
		M:createTextures(0, "M", 1), //Meutriere
		B:createTextures(0, "windows", 7, 240,20),
		R:createTextures(0, "windows", 8),
		x:createTextures(0, "windows", 10), //Rock
		w:[
			createTextures(0, "windows", 1),
			createTextures(0, "windows", 2),
			createTextures(0, "windows", 3),
			createTextures(0, "windows", 4),
			createTextures(0, "windows", 5),
			createTextures(0, "windows", 6)
		],
		b:[]
	}
	for (var i=0; i<=50; i++) {
		textures.b[i] = createTextures(i);
	}
	
	//
	// Init pseudo random for mountain and decorations
	//
	mr = new random(201);
	generateMountain();
	generateMountain();
	generateMountain();
};
var goHome;
function drawPlane(){
	ctx.beginPath();
	ctx.strokeStyle = 'black'; 
	ctx.moveTo(-w, 0.5*zoom);
	ctx.lineTo( w, 0.5*zoom);
	ctx.stroke();
};
function fillForm(ax,ay,aw,ah, radius) {
	ctx.beginPath();
	if (radius) {
		ctx.arc(ax*zoom, -ay*zoom, radius/2 * zoom, 0, 2 * Math.PI);
	} else {
		ctx.rect(ax*zoom, -ay*zoom, aw*zoom, ah*zoom);
	}
	
	ctx.translate(ax*zoom, ay*-zoom);
	ctx.scale(zoom/40, zoom/40);
	ctx.fill();
	ctx.stroke();
	ctx.scale(1/(zoom/40), 1/(zoom/40));
	ctx.translate(ax*-zoom, ay*zoom);
}
var font = [];
font[1] = "010110010010111";
font[2] = "111001111100111";
font[3] = "111001111001111";
font[4] = "101101111001001";
font[5] = "111100111001111";
font[6] = "111100111101111";
font[7] = "111001010010010";
font[8] = "111101111101111";
font[9] = "111101111001111";
font["A"] = "111101111101101";
font["B"] = "110101111101111";
font["C"] = "111100100100111";
font["D"] = "110101101101110";
font["E"] = "111100111100111";
font["I"] = "111010010010111";
font["K"] = "101101110101101";
font["L"] = "100100100100111";
font["M"] = "101111111101101";
font["O"] = "111101101101111";
font["P"] = "110101110100100";
font["R"] = "111101111110101";
font["S"] = "111100111001111";
font["T"] = "111010010010010";
font["U"] = "101101101101111";
font["V"] = "101101101101010";

function drawText(posX, posY, text, fontSize=10, fillStyle) {
	ctx.beginPath();
	ctx.fillStyle = fillStyle || textures.T;
	pX = posX;
	for (var i=0; i<text.length; i++) {
		var index = 0;
		for (var y=0; y<5; y++) {
			for (var x=0; x<3; x++) {
				if (text[i]==" ") { continue; }
				if (font[text[i]][index++]=="1") {
					ctx.rect(
						pX+x*fontSize,
						posY+y*fontSize,
						fontSize,
						fontSize
					);
				}
			}
		}
		pX += fontSize*3.5; //Interlettrage
	}
	ctx.save();
	ctx.translate(posX,posY);
	ctx.scale(0.5,0.5);
	ctx.fill();
	ctx.restore();
}
function render(){
	// Clear the canvas
	
	ctx.fillStyle = "#4290D1"; //Sky color
	ctx.fillRect(0,0,w,h);
	
	//
	// Move clouds
	//
	for (i=0; i<clouds.length; i++) {
		if (clouds[i].x<-200) { clouds[i].x = 1200; }
		ctx.drawImage(cloud, clouds[i].x-=clouds[i].speed,clouds[i].y, 128*clouds[i].size, 71*clouds[i].size);
	}
	
	//
	// Les montagnes
	//
	ctx.drawImage(mountain,0,0,w,h);
	
	if (screen=="home") {
		drawText(372,75,"CASTLE ATTACK");

		x=180;
		for (i=0; i<7; i++) {
			ctx.beginPath();
			ctx.fillStyle = "white";
			ctx.rect(x, 250, 80, 60);
			ctx.strokeStyle = i>localStorage.gp44.length ? "gray" : "black";
			ctx.fill();
			ctx.stroke();
			
			ctx.beginPath();
			drawText(x+18,260,"LEVEL "+(i+1), 2., i>localStorage.gp44.length ? "gray" : "black");

			for (var j=0; j<3; j++) {
				ctx.globalAlpha = localStorage.gp44[i]>j ? 1 : .2;
				ctx.drawImage(images["windows9"], x+12+(20*j),290,16,12);
			}
			ctx.globalAlpha = 1;
			x+=130;
		}
		ctx.beginPath();
	}
	if (screen=="level") {
		//
		// Les ammos
		//
		if (!goHome) {
			var m = areas[levelNumber].ammos; //Toutes les ammos du niveau
			var inProgress = m.length-ammos.length;

			var y = 20;
			for (l=0; l<m.length; l++) {
				drawText(30, y-(l==inProgress ? 4 : 0),m[l], l==inProgress ? 3 : 2, l<inProgress ? "#666" : null);
				y+=20;
			}

			ctx.beginPath();
			ctx.strokeStyle = "black";
			if (ammos.length) {
				ctx.moveTo(10, 22+(inProgress*20));
				ctx.lineTo(25, 22+(inProgress*20));
			} else {
				end();
			}
			ctx.stroke();
		}
		
		//
		// Les éléments 
		//
		ctx.save();
		ctx.translate(w/2 + camera.x, h/2 + camera.y);  // Translate to the center
		drawAll();
		drawPlane();
		ctx.restore();

		//
		// Le Baliste
		//
		ctx.fillStyle = textures.B;
		ctx.save();
		ctx.translate(w/2 + camera.x - 15*zoom, h/2 + camera.y - 0*-zoom);

		// La barre de lancement
			ctx.translate(0, 0.25*-zoom);
			ctx.rotate(-weaponAngle);
				fillForm(0,0.25,4,0.5);

				// Viseur
				ctx.fillRect(3.5*zoom, -0.*zoom, 0.25*zoom, -0.5*zoom);

				// Le fil
				var puissance = keys[32] * ((Math.min(frames-power, 100))/100);

				ctx.beginPath();
				ctx.moveTo(3.5*zoom, -0.5*zoom);
				ctx.lineTo((3.5-3.5*(puissance))*zoom, -0.25*zoom);
				ctx.stroke();
			ctx.rotate(weaponAngle);
			ctx.translate(0, 0.25*zoom);

		// Le socle
			fillForm(-0.5,0.25,6,0.5);

		// Les roues
			fillForm(0.5,0,1,1,1);
			fillForm(4.5,0,1,1,1);

		ctx.restore();
	}
	
};

screen = "";
function gotoScreen(s,i) {
	console.log("gotoScreen", s);
	screen = s;
	if (screen=="level") {
		levelNumber=i;
		initLevel(i); 
	}
}
gotoScreen("home");


var weaponAngle = 0.1;
var frames = 0;
var explode = {frame:0, position:{}};
var camera = {x:-160, y:190};
var zoom = 20;
// Animation loop
function animate(){
	if (explode.frame) {
		explode.frame++;
		explodeForce = (1*(30-explode.frame));
		camera.x=explode.position.x + (Math.random() * explodeForce)-(explodeForce/2);
		camera.y=explode.position.y + (Math.random() * explodeForce)-(explodeForce/2);
		if (explode.frame==30) {
			explode.frame = 0;
			camera.x = explode.position.x;
			camera.y = explode.position.y;
		}
	}
	for (var i=0; i<balls.length; i++) {
		var ball = balls[i];
		if (ball.weapon== "TIME BOMB") {
			if (frames - ball.age>300) {
				ball.explode();
			}
		}
	}
	frames++;
	if (keys[38]) {
		weaponAngle = Math.min(1.17, weaponAngle+0.015);
	}
	if (keys[40]) {
		weaponAngle = Math.max(0.1, weaponAngle-0.015);
	}
	if (keys[27]) {
		gotoScreen("home");
	}
	requestAnimationFrame(animate);

	// Move physics bodies forward in time
	if (screen=="level" && world) {
		world.step(1/60);
	}

	// Render scene
	render();
};


onclick = function(e) {
	if (screen=="home") {
		var rect = canvas.getBoundingClientRect();

		// Récupérez la position du clic par rapport au coin supérieur gauche du canvas (non redimensionné)
		var x = e.clientX * canvas.width / rect.width;
		var y = e.clientY * canvas.height / rect.height;

		xTest=180; //180->260

		if (y>250 && y<310) {
			for (i=0; i<7; i++) {
				if (x>=xTest && x<=xTest+80) {
					if (i<=localStorage.gp44.length) {
						gotoScreen("level",i)
					}
				}
				xTest+=130;
			}
		}
	}
}
onkeyup = function(e) {
	keys[e.keyCode] = false;
	if (e.keyCode==32) {
		if (ammos.length) {
			var force = Math.sqrt((Math.min(frames-power, 100))/100*20)*4;
			power = 0;

			//Balle
			var balle = addObject({
				type:"circ",
				weapon:ammos.shift(),
				radius:0.20,
				mass:1,
				position:[
					-15 + Math.cos(weaponAngle)*3.5,
					0.25 + Math.sin(weaponAngle)*3.5
				]
			});
			world.gravity = vec2.fromValues(0, -5); //-0.978);
			balle.age = frames;
			balle.body.vel = [Math.cos(weaponAngle)*force,Math.sin(weaponAngle)*force];
			balls.push(balle);
		}
	}
};

var power = 0;
onkeydown = function(e) {
	keys[e.keyCode] = true;
	if (!power && e.keyCode==32) {
		power = frames;
	}
}

localStorage.gp44 = localStorage.gp44 || "";
init();
animate();

function end() {
	if (!goHome) {
		goHome = true;
		setTimeout(function() {
			localStorage.gp44 = localStorage.gp44.slice(0, levelNumber) + (Math.max(localStorage.gp44[levelNumber]|0,3-nbChests)) + localStorage.gp44.slice(levelNumber + 1);
			
			gotoScreen("home");
			goHome = false;
		}, 10000);
	}
}