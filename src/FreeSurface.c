/*
 * FreeSurface.c
 *
 *  Created on: Jan 20, 2017
 *      Author: abauville
 */


void FreeSurface_init() {

}

void FreeSurface_advectErodeDepose() {
	// Pseudo code for the to advect and redefine the marker chain:
	/*

	//Advect
	for i {
		xTemp = FreeSurface->x[i];
		yTemp = FreeSurface->y[i];
		vx = stuff;
		vy = stuff;

		xTemp += vx*dt;
		yTemp += vy*dt;

	}

	// Reinit the y positions
	for all nodes of the marker chain{
		y = 0.0;
	}
	// The interp1D  for the fixed x (on nodes or cell centers)
	for segment {
		segment_x0 = stuff;
		segment_x1 = stuff;
		segment_y0 = stuff;
		segment_y1 = stuff;
		if segment_x1>segment_x0 {
			// find which node of the marker chain is contained in between the two ends of the segment
			for (nodes inside) {
				xTemp = linearinterp on the segment;
				yTemp = linearinterp on the segment;
				if (yTemp>y) {
					y = yTemp;
					x = xTemp;
				}
			}
		} else {
			// means there is a zigzag, there must be another segment above that follows segment_x1>segment_x0
			// do nothing
		}
	}



	// Apply erosion and sedimentation if needed
	for 1:end-1 {
		// 1d diffusion of the y position of marker

		// sedimentation by changing the y position of markers
	}

	// Change particles phase
	for relevant nodes {
		// check particles if above or below the marker chain
		// if air or water < surface
		// change to something

		// if solid > surface
		// change to air or water


	}


*/

}
