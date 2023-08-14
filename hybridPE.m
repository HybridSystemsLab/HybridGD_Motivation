classdef hybridPE < HybridSystem
    
    properties(SetAccess=immutable)
        theta
        gammac
        gammad
        lambdac
        lambdad
    end

    properties 
        nu
        R
    end
    
    methods 
        function this = hybridPE(parameters)
            state_dim = 11; % (x,thetahat,psi,eta,u)
            this = this@HybridSystem(state_dim);
            
            this.theta = parameters.theta;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
            this.lambdac = parameters.lambdac;
            this.lambdad = parameters.lambdad;
            this.nu = parameters.nu;
        end

        function zdot = flowMap(this, z, t, j)            
            x = z(1:2);
            thetahat = z(3:4);
            psi = reshape(z(5:8),2,2);
            eta = z(9:10);
            u = z(11);
            
            m = this.nu*sin(2*t)*[1; 1];
            xm = x + m;

            f_c = [0; 0];
            phi_c = [sin(t) 0; 0 0];
            xdot = f_c + phi_c*this.theta;

            if this.R == 1
                thetahatdot = 0*thetahat;
                etadot = 0*eta;
                psidot = 0*psi;
            else
                psidot = -this.lambdac*psi + phi_c;
                etadot = -this.lambdac*(xm + eta) - f_c;
                y = xm + eta;
                thetahatdot = this.gammac*psi.'*(y - psi*thetahat);
            end

            udot = 1;

            zdot = [xdot; thetahatdot; psidot(:); etadot; udot];
        end

        function zplus = jumpMap(this, z, t, j)
            x = z(1:2);
            thetahat = z(3:4);
            psi = reshape(z(5:8),2,2);
            eta = z(9:10);
            u = z(11);
            
            g_d = [0; 0];
            % phi_d = (pi/3)*[1 2; 2 4];
            phi_d = [1 2; 2 4];
            xplus = g_d + phi_d*this.theta;

            m = this.nu*sin(2*t);
            xm = x + m;
            xmplus = xplus + m;

            if this.R == 0
                thetahatplus = thetahat;
                etaplus = eta;
                psiplus = psi;
            else
%                 y_d = xmplus - f_d;
%                 thetahatplus = thetahat + (phi_d.'/(this.gammad + norm(phi_d)^2))*(y_d - phi_d*thetahat);
%                 psiplus = (1 - this.lambdad)*psi;
%                 etaplus = (1 - this.lambdad)*(xm + eta) - xmplus;
                
                psiplus = (1-this.lambdad)*psi + phi_d;
                etaplus = (1-this.lambdad)*(xm + eta) - g_d;
                yplus = xmplus + etaplus;
                thetahatplus = thetahat + (psiplus.'/(this.gammad + norm(psiplus)^2))*(yplus - psiplus*thetahat);
            end

            uplus = 0;

            zplus = [xplus; thetahatplus; psiplus(:); etaplus; uplus];
        end
        
        function inC = flowSetIndicator(this, z, t, j)
            x = z(1:2);
            thetahat = z(3:4);
            psi = reshape(z(5:8),2,2);
            eta = z(9:10);
            u = z(11);
            
            if u <= 2*pi
                inC = 1;
            else
                inC = 0;
            end
        end

        function inD = jumpSetIndicator(this, z, t, j)
            x = z(1:2);
            thetahat = z(3:4);
            psi = reshape(z(5:8),2,2);
            eta = z(9:10);
            u = z(11);
            
            if u >= 2*pi
                inD = 1;
            else
                inD = 0;
            end
        end
    end
end
