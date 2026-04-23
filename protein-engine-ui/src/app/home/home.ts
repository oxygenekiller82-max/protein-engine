import { Component,ChangeDetectorRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Router,RouterLink } from '@angular/router';
import { HttpClient, HttpClientModule } from '@angular/common/http';

import { SequenceService } from '../services/sequence.service';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [CommonModule, FormsModule,RouterLink,HttpClientModule],
  templateUrl: './home.html',
  styleUrl: './home.css'
  
})
export class HomeComponent {
  

  //FOR TESTING REMOVE THESE TODO 
  hydrophobicityMin: number =0.0;
  hydrophobicityMax: number =18.0
  massMin: number =10500;
  massMax: number =11500;
  stabilityMin: number =80.0;
  bindingMin: number =220.0;

  sequenceName: string = '';
  animationText: string = '';
  generated: boolean = false;

  //hydrophobicityMin: number | null = null;
  //hydrophobicityMax: number | null = null;
  //massMin: number | null = null;
  //massMax: number | null = null;
  //stabilityMin: number | null = null;
  //bindingMin: number | null = null;

  //sequenceName: string = '';
  //animationText: string = '';
  //generated: boolean = false;

  //PLACE HOLDER
  public finalSequence: string[] = [];


  //flaks screams 
  public errorMessage: string = '';
  

  //XAI YES
  //behaviour subject my goat 
  history: any[] = [];
  currentFrame: number = 0;
  displayedSequence: string[] = [];
  lastAction: string = '';
  lastPruneReason: string = '';
  isAnimating: boolean = false;

   //Analyze seqeunce thingy 
   peptideStats: any= null;
    

  //new for the biological switch and the target length: 

  length: number= 20; //default 20 i guess ??
  biologicalSwitch: boolean = true;

   constructor(private router: Router,
    private seqService: SequenceService,
    //this wakes up Agnualr jsut in case the UI doesn't auto update WHATTT
    private cdr: ChangeDetectorRef,
    private http: HttpClient 
  ) {}

  //XAI real deal YUH
  startAnimation(history: any[]) {
    this.history = history;
    this.currentFrame = 0;
    this.isAnimating = true;
    this.playNextFrame();
  }

  playNextFrame(){
    if(this.currentFrame>=this.history.length){
      this.isAnimating=false;
      return;  //WHAT ???? TODO
    }
    //else 
    const frame=this.history[this.currentFrame];
    this.lastAction=frame.action;
    this.lastPruneReason=frame.prune_reason;
    this.displayedSequence=[...frame.sequence]; //WHAT ? TODO
    this.currentFrame++;
    this.cdr.detectChanges();

    setTimeout(() => this.playNextFrame(), 150);//150ms 
    //this should have a button i think speed up or down!! TODO

  }

  
   //really ? c:
   downloadLog(): void{
    if(!this.history.length) return; 

    const header = 'STEP     | ACTION     | LVL    | LAST_AA    | PRUNE_REASON         | SEQUENCE\n' +
                 '--------------------------------------------------------------------------------------------------------------\n';

    const rows= this.history.map(frame=>{
      const step = String(frame.step).padEnd(8);
      const action = String(frame.action).padEnd(11);
      const level = String(frame.level).padEnd(7);
      const lastAA = String(frame.last_AA_added).padEnd(11);
      const reason = String(frame.prune_reason).padEnd(36);
      const sequence = frame.sequence.join('-');

      return `${step} | ${action} | ${level} | ${lastAA} | ${reason} | ${sequence}`;
    }).join('\n');

    //blob API -> creates file from memory + triggers download 
    const blob = new Blob([header+rows],{type: 'text/plain'});
    const url=URL.createObjectURL(blob);
    const a = document.createElement('a');

    a.href=url;
    a.download='SEARCH_LOG.log';
    a.click();
    URL.revokeObjectURL(url);
  }
  


   //the real deal YUH
  generate(): void {
    //alright no more Nulls. 
    const targets = [
      this.hydrophobicityMin, this.hydrophobicityMax, 
      this.massMin, this.massMax, 
      this.stabilityMin, this.bindingMin, 
      this.length
    ];

    //.some lool 
    if (targets.some(value => value === null || value === undefined)) {
      this.generated = true;
      this.animationText = "Validation Error";
      this.errorMessage = "Please fill in all biological constraints before generating.";
      this.cdr.detectChanges();
      return; // STOP HERE - Don't call the service!
      }

    //else all filled 
    this.generated = true;
    this.animationText = 'Generating...';
    this.errorMessage = '';

    //loading state ?
    this.generated = true;
    this.animationText = 'Generating the suitable sequence...';

    const payload={
      target_length: this.length,
      biological_switch: this.biologicalSwitch,
      targets: {
        hydro_min: this.hydrophobicityMin,
        hydro_max: this.hydrophobicityMax,
        mass_min: this.massMin,
        mass_max: this.massMax,
        stability_min: this.stabilityMin,
        binding_min: this.bindingMin,
        }
    };

    //service call with that payload YES
    this.seqService.generate(payload);
    

    //once result arrives -> change text 
    this.seqService.currentResult$.subscribe(data=>{
      if(data){

        if(data.status !='done'){
          this.generated = true;
          this.animationText = "Error:";
          this.errorMessage = data.message;
          this.finalSequence = [];
        }else{
          this.errorMessage="";
          this.finalSequence = data.sequence;
          this.startAnimation(data.history);
          //NOW the stats data.sequence is the input 

          this.http.post<any>('http://localhost:8081/api/PFA/sequence_stats', {
            
            sequence: data.sequence
                }).subscribe({
                    next: (stats) => {
                        this.peptideStats = stats;
                        this.cdr.detectChanges();
                    },
                    error: (err) => console.error('Stats failed', err)
                });
        }

        this.generated=true;
        this.animationText="DONE GENERATING!!!"
        //TEMPRORAY i wanna see it 
        console.log("LOG:",data);
        //ALSO TEMP saving to arary jsut to see it 
        //this.finalSequence = data.sequence;
        //TEMPORARY TODO 
        

        this.cdr.detectChanges();// ANGULAR WAKE THE FICK UP !! 
        //REDIRECT TODO 
      }
    });

   // setTimeout(() => {
    //  this.animationText = 'Séquence générée avec succès ';
    //}, 1500);

  }

  analyser(): void {
    if (!this.generated) {
      alert('Veuillez d\'abord générer une séquence.');
      return;
    }
    alert('Analyse lancée !');
  }

    compare(): void {
    this.router.navigate(['/compare']);
  }
}