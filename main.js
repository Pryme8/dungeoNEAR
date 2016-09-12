// Dungeonear v0.0.1
// Experments in Procedural Level Creation
// Author: Andrew V Butt Sr.
// Pryme8@gmail.com || http://pryme8.github.io

dungeon = function(args){
	args = args || {};
	args.zSize = args.zSize || 60;
	args.zDiv = args.zDiv || 3;
	if(typeof args.baseNoise === 'undefined' || typeof args.baseNoise !== 'object' ){
	return "Error need a base Noise from Das Noise!";
	}
	
	this.noise = args.baseNoise;
	this.zSize= args.zSize;
	this.zDiv = args.zDiv;
	this._createRefMap();
};

dungeon.prototype._createRefMap = function(){
	var cells = [];
	var pCount = this.zDiv * this.zDiv;
	
	function perm(s,c){
		if (c == 0) {
        cells.push(s);
        return;
    	}
		perm(s+'0', c-1);
		perm(s+'1', c-1);
	}
	
	perm('',pCount);
	
	var last = cells.splice(cells.length-1,1);
	cells.unshift(last+'');
	//console.log(cells);
	
	
	var map = [];
	
	for(var i = 0; i < cells.length; i++){
		map.push(new dungeon.Zone(this.zSize, this.zDiv, cells[i]));	
	}
	this.map = map;
	
	//console.log(map);
	var self = this;
	setTimeout(function(){self._calculateMap();},0);
	
};

dungeon.prototype._calculateMap = function(){
	var map = this.map;
	var cvas = document.createElement('canvas');
	var ctx = cvas.getContext('2d');
	
	var X = 0, Y = 0;
	var cellSize = this.zSize/this.zDiv;
	
	cvas.width = 20*this.zSize;
	cvas.height = Math.ceil(map.length/20)*this.zSize;
	
	for(var i = 0; i < map.length; i++){
		var x = 0, y = 0;
		for(var j = 0; j < map[i].cells.length; j++){
			
			if(map[i].cells[j] == 1){
				ctx.fillStyle = "#FFF";
			}else{
				ctx.fillStyle = "#000";
			}
			
			ctx.fillRect(x+X,y+Y,cellSize,cellSize);

			x+=cellSize;
			if(x > this.zSize-cellSize){
				y+=cellSize;
				x=0;
			}
		};
		
		ctx.strokeStyle = "rgba(255,0,0,0.2)";	
		ctx.strokeRect(X,Y,this.zSize,this.zSize);
		
		var imgData = ctx.getImageData(X,Y,this.zSize,this.zSize);
		
		map[i].imgData = imgData;
		map[i].x = X;
		map[i].y = Y;
				
		X+=this.zSize;
		if(X > cvas.width-this.zSize){
			Y+=this.zSize;
			X=0;
		}
	};	
};

dungeon.prototype._idNoise = function(x,y,noise){
	if(typeof noise ==='undefined'){
	noise = this.noise;	
	}
	var cellSize = this.zSize/this.zDiv;
	
	var string = '';
	var self = this;
	var ctx = (document.getElementById('noise-canvas')).getContext('2d');
		ctx.fillStyle = "red";
		
		var cX = 0;
		var cY = 0;
	
	for(var i=0; i<this.zDiv*this.zDiv; i++){
		var t = 0;

		for(var pY = 0; pY < cellSize; pY++){
			for(var pX = 0; pX < cellSize; pX++){
				t+=noise.getValue({x:(pX+(this.zSize*x)+(cellSize*cX)),
								   y:(pY+(this.zSize*y)+(cellSize*cY))});
									
			}
		}
		
		t/=(cellSize*cellSize);
		if(t<0.45){
		string+=0+"";
		}else{
		string+=1+"";
		}

		cX++;
		if(cX>this.zDiv-1){
		cX=0;
		cY++;
		}
	}

	for(var i=0; i<this.map.length; i++){
		if(this.map[i].searchString == string){
			return this.map[i];
		}
	};
};

dungeon.Zone = function(size, div, state){
	this.size = size;
	this.div = div;
	this.searchString = state;
	this.cells = state.split('');
	for(var i = 0; i<this.cells.length; i++){
	this.cells[i] = parseFloat(this.cells[i], 10);	
	}
	return this;	
};