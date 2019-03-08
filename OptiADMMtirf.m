classdef OptiADMMtirf < OptiADMM
    % OptiADMMtirf: Derivation of the OptiADMM class that includes a step
    % for background estimation for MA-TIRF reconstruction as proposed in [1].
    %
    % :param F_0: see OptiADMM
    % :param F_n: see OptiADMM
    % :param H_n: see OptiADMM
    % :param rho_n: see OptiADMM
    % :param estiback: Bollean (true to activate background estimation)
    % :param lambBack: Regularization parameter for the background estimation
    %
    % **Reference**
    %
    % [1] Emmanuel Soubies, Laure Blanc-Feraud, Sebastien Schaub and Ellen Van Obberghen-Schilling,
    %     "Improving 3D MA-TIRF Reconstruction with Deconvolution and
    %     Background Estimation", Proc ISBI 2019
    %
    % **Example** ADMM=OptiADMM(F0,Fn,Hn,rho_n,solver)
    %
    % See also :class:`Opti`, :class:`OptiADMM`, :class:`OutputOpti`, :class:`Cost`
    
    %%    Copyright (C) 2018
    %     E. Soubies emmanuel.soubies@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    

    % Full public properties
    properties
        yy;
        Op;
        back;
        estiback
        sz
    end
    
    methods
        %% Constructor
        function this=OptiADMMtirf(F0,Fn,Hn,rho_n,estiback,lambBack)
            if nargin <6 || isempty(lambBack), lambBack=0; end;
            this@OptiADMM(F0,Fn,Hn,rho_n);
            this.name='OptiADMMtirf';
            this.yy=Fn{1}.H1.y;
            this.sz=Fn{1}.H1.sizein;
            G=LinOpGrad(this.sz(1:2));
            this.Op=LinOpIdentity(this.sz(1:2))+lambBack*(G'*G);
            this.estiback=estiback;
        end
        function updateParams(this)
            % Update Background signal
            if this.estiback 
                this.back=max(this.Op.applyInverse(sum(this.yy-this.Fn{1}.H2*this.yn{1},3)/sqrt(this.sz(3))),0);
                this.Fn{1}.H1.y=this.yy-repmat(this.back,[1,1,size(this.yy,3)])/sqrt(this.sz(3));
            end
        end
    end
end
