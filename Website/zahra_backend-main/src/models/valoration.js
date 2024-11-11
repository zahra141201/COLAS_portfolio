'use strict';
const { Model } = require('sequelize');

module.exports = (sequelize, DataTypes) => {
  class Valoration extends Model {
    static associate(models) {
      this.belongsTo(models.User, {
        foreignKey: 'email_user',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les évaluations associées lorsque l'utilisateur évalué est supprimé
      });
      this.belongsTo(models.User, {
        foreignKey: 'made_by',
        onDelete: 'CASCADE',
        onUpdate: 'CASCADE', // Supprime les évaluations associées lorsque l'utilisateur qui a fait l'évaluation est supprimé
      });
    }
  }

  Valoration.init({
    comment: DataTypes.STRING,
    puntuation: DataTypes.INTEGER,
    email_user: DataTypes.STRING,
    made_by: DataTypes.STRING
  }, {
    sequelize,
    modelName: 'Valoration',
  });

  return Valoration;
};
