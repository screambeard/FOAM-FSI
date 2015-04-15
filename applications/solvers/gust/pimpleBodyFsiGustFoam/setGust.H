const scalar& time = runTime.time().value();

Ugust == gustObj.getGustVelocity(time,U);
phiGust == (Ugust & mesh.Sf());