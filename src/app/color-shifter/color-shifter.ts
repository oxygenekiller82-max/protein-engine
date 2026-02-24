import { Component } from '@angular/core';

@Component({
  selector: 'app-color-shifter',
  imports: [],
  templateUrl: './color-shifter.html',
  styleUrl: './color-shifter.css',
})
export class ColorShifter {
  //DEFINE STATE !!
  bgColor: string ='#f0f0f0';
  clickCount: number =0;

  //LOGIC HERE
  changeColor(){
    const colors=['#ff5733', '#33ff57', '#3357ff', '#f333ff', '#33fff3'];
    const randomIndex = Math.floor(Math.random()*colors.length);

    this.bgColor = colors[randomIndex];
    this.clickCount++;

    console.log(`Color changed, new color =${this.bgColor}. Clicked ${this.clickCount} times.`);
  }

}
